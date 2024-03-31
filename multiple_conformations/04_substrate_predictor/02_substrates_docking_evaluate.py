#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import logging
import json
import functools

from libs import utils
from libs import paths
from libs import data_model
from libs.configuration import initialize_configuration, TASK, ENV
from libs.configuration import _SPECIAL_RESIDUES
from libs import docking_evaluate
from libs.data_model import SubstrateEntry, Fragment, MappingEntry, McsaEntry
from libs import file_formats
from libs import active_site as active_site_service

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration()
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        substrates_docking_evaluate()
        logger.info("Execution time: %d s", time.time() - start_time)


def substrates_docking_evaluate():
    logger.info("Evaluating docking '%s' -> '%s'", paths.complexes(),
                paths.substrates_eval())
    task_info = _read_task_info()
    _check_task_info(task_info)
    database = data_model.create_read_only_from_directory(
        paths.knowledge_base_path())
    in_dir = paths.complexes()
    out_path = paths.substrates_eval()
    number_of_fragments = 0
    with open(out_path, "w") as out_stream:
        dockfiles = os.listdir(in_dir)
        dockfiles.sort()
        for file_name in dockfiles:
            if not file_name.endswith(".pdbqt"):
                continue
            name = file_name[:file_name.rfind(".")]
            try:
                fragments = _evaluate_file(
                    database, task_info["fragment"], in_dir, name)
            except RuntimeError:
                logger.exception("Can't evaluate file: %s", name)
                continue
            _remove_duplicate_atom_pairs(fragments)
            number_of_fragments += len(fragments)
            data_model.dump_entry(fragments, out_stream)
            out_stream.write("\n")
        logger.info("Output %d fragments.", number_of_fragments)


def _read_task_info():
    with open(paths.task_info()) as in_stream:
        return json.load(in_stream)


def _check_task_info(task_info):
    if "fragment" not in task_info:
        message = "Please insert fragment information in task.json file."
        raise RuntimeError(message)


def _evaluate_file(database, fragment_definition, directory, name):
    containers = _load_fragments(database, fragment_definition, directory, name)
    start_count = len(containers)
    docking_evaluate.evaluate_fragments_basic(database, containers)
    containers = docking_evaluate.filter_out_bad(containers)
    docking_evaluate.evaluate_fragments_expensive(database, containers)
    logger.debug("Loaded %d (%d before filtering) fragments from '%s'", len(containers), start_count, name)
    return [container.fragment for container in containers]


def _load_fragments(database, fragment_definition, directory, name):
    info_path = os.path.join(directory, name + ".json")
    with open(info_path) as in_stream:
        info = json.load(in_stream)

    substrate = database.read(SubstrateEntry, info["ref"]["substrate"])
    pdbqt_path = paths.substrate_pdbqt_path(substrate)
    sdf_path = paths.substrate_sdf_path(substrate)

    docking_path = os.path.join(directory, name + ".pdbqt")

    docking_result = docking_evaluate.load_docking_result(
        pdbqt_path, sdf_path, docking_path)

    fragments = []
    for model_index, molecule in enumerate(docking_result):
        discovered_fragment = docking_evaluate.extract_fragments(
            fragment_definition, molecule)
        for frag_index, atoms in enumerate(discovered_fragment):
            fragment = _create_fragment(
                info, fragment_definition, model_index,
                frag_index, molecule, atoms)
            fragments.append(docking_evaluate.FragmentContainer(
                fragment, molecule, atoms, fragment_definition["names"]))
    return fragments


def _create_fragment(
        info, fragment_definition, model_index, frag_index, molecule, atoms):
    fragment = Fragment()
    fragment.mapping = info["ref"]["mapping"]
    fragment.substrate = info["ref"]["substrate"]
    fragment.name = info["name"]
    fragment.model_index = model_index
    fragment.fragment_index = frag_index
    fragment.vina_score = molecule.prop("vina")
    fragment.atoms = []
    for atom, name in zip(atoms, fragment_definition["names"]):
        fragment_atom = Fragment.Atom()
        fragment_atom.x = atom.x
        fragment_atom.y = atom.y
        fragment_atom.z = atom.z
        fragment_atom.serial = atom.serial
        fragment_atom.name = name
        fragment.atoms.append(fragment_atom)
    return fragment


@functools.lru_cache(maxsize=16)
def load_protein(database, mapping_ref):
    """
    :return Atoms that are part of active site and are specified in MCSA record.
    """
    mapping = database.read(MappingEntry, mapping_ref)
    mcsa = database.read(McsaEntry, mapping.mcsa().first().value())

    protein_path = paths.protein_pdbqt_path(mapping, mcsa)
    protein_pdbqt = file_formats.load_pdbqt_file(protein_path)

    active_site = active_site_service.active_site_from_mcsa(mcsa)
    atoms_by_names = active_site_service.reacting_atoms(
        protein_pdbqt, mapping, active_site)

    reactive_atoms = [atom for atoms in atoms_by_names.values()
                      for atom in atoms]

    return {
        "reactive_by_name": atoms_by_names,
        "reactive": reactive_atoms
    }


# noinspection SpellCheckingInspection
def _remove_duplicate_atom_pairs(fragments):
    for fragment in fragments:
        to_evaluate = []
        selected_distances = []
        for distance in fragment.distances:
            if (distance["name"][1] in _SPECIAL_RESIDUES) and (distance["name"][3] != "NZ"):
                to_evaluate.append(distance)
            else:
                selected_distances.append(distance)
        if len(to_evaluate) > 0:
            fragment.distances = selected_distances
            evaluated_distances = set()
            for i in range(len(to_evaluate)):
                best_distance = to_evaluate[i]
                resi, resn, role, atom, atype = best_distance["name"]
                if (resi, resn, role, _SPECIAL_RESIDUES[resn][atom], atype) in evaluated_distances:
                    continue
                _dist = {"name": None, "val": None, "old_name": None}
                for j in range(i+1, len(to_evaluate)):
                    _resi, _resn, _role, _atom, _atype = to_evaluate[j]["name"]
                    if resi==_resi and resn==_resn and role==_role and atype==_atype:
                        _dist["name"] = (resi, resn, role, _SPECIAL_RESIDUES[resn][atom], atype)
                        if best_distance["val"] < to_evaluate[j]["val"]:
                            _dist["val"] = best_distance["val"]
                            _dist["old_name"] = atom
                        else:
                            _dist["val"] = to_evaluate[j]["val"]
                            _dist["old_name"] = _atom
                            best_distance = to_evaluate[j]
                fragment.distances.append(_dist)
                evaluated_distances.add((resi, resn, role, _SPECIAL_RESIDUES[resn][atom], atype))


docking_evaluate.load_protein = load_protein

if __name__ == "__main__":
    main()
