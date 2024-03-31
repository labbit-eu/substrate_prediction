#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import logging
import json
import functools
import multiprocessing as mp
import gzip

from libs import utils
from libs import paths
from libs import data_model
from libs.configuration import initialize_configuration, TASK, ENV
from libs.configuration import _SPECIAL_RESIDUES
from libs import docking_evaluate
from libs.data_model import Fragment, MappingEntry, McsaEntry
from libs import file_formats
from libs import active_site as active_site_service

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration()
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    conformers = 1
    with utils.FileLogging(paths.logs(), __file__):
        for c in range(conformers):
            if utils.among_candidate_processing_statuses(["docking-failed", "docking-completed"], c):
                candidates_docking_evaluate(c)
            else:
                logger.error("No results to evaluate, run EXT09_candidates_docking.py first")
            logger.info("Execution time: %d s", time.time() - start_time)


def candidates_docking_evaluate(conformer):
    logger.info("Evaluating docking '%s' in %d parallel processes", paths.candidates(), ENV["num_parallel_cpu"])
    candidates_info = _read_candidates_info(conformer)
    task_info = _read_task_info()
    _check_task_info(task_info)
    database = data_model.create_read_only_from_directory(
        paths.knowledge_base_path())
    mapping_ref = candidates_info["mapping"]

    pool = mp.Pool(processes=ENV["num_parallel_cpu"])
    processing = []
    for index, files in enumerate(utils.split_files_in_subdirs_to_sets(paths.candidates_conf_docking(conformer))):
        processing.append(pool.apply_async(_candidates_docking_evaluate_sets, args=(index, files, database, task_info,
                                                                                    mapping_ref, conformer)))

    numbers = [0, 0]

    for p in processing:
        a = p.get()
        numbers[0] += a[0]
        numbers[1] += a[1]
    pool.close()

    logger.info("Output %d fragments extracted from %d files.", numbers[0], numbers[1])


def _candidates_docking_evaluate_sets(index, files, database, task_info, mapping_ref, conformer):
    out_path = paths.candidates_eval(index)
    if os.path.isfile(out_path):
        writing_mode = "at"
    else:
        writing_mode = "wt"
    number_of_fragments = 0
    num_files = 0
    if TASK["candidates"]["docking_mode"] == "manual":
        files = _merge_dockings(files)
    with gzip.open(out_path, writing_mode) as out_stream:
        for docking_file in files:
            if not docking_file.endswith(".pdbqt"):
                continue
            try:
                fragments = _evaluate_file(database, task_info["fragment"], mapping_ref, docking_file, conformer)
            except RuntimeError:
                logger.exception("Can't evaluate file: %s", docking_file)
                continue
            _remove_duplicate_atom_pairs(fragments)
            number_of_fragments += len(fragments)
            num_files += 1

            data_model.dump_entry(fragments, out_stream)
            out_stream.write("\n")

    logger.info("%d files were evaluated.", num_files)
    return [number_of_fragments, num_files]


def _read_candidates_info(conformer):
    with open(paths.candidates_conf_info(conformer)) as in_stream:
        return json.load(in_stream)


def _read_task_info():
    with open(paths.task_info()) as in_stream:
        return json.load(in_stream)


def _check_task_info(task_info):
    if "fragment" not in task_info:
        message = "Please insert fragment information in task.json file."
        raise RuntimeError(message)


def _evaluate_file(database, fragment_definition, mapping_ref, docking_file, conformer):
    containers = _load_fragments(fragment_definition, mapping_ref, docking_file)

    initial_fragments_count = len(containers)
    docking_evaluate.evaluate_fragments_basic(database, containers, conformer)
    containers = docking_evaluate.filter_out_bad(containers)
    docking_evaluate.evaluate_fragments_expensive(database, containers, conformer)

    logger.debug("Loaded %d (%d before filtering) fragments from '%s'", len(containers), initial_fragments_count,
                 docking_file)
    return [container.fragment for container in containers]


def _load_fragments(fragment_definition, mapping_ref, docking_file):
    file_name = os.path.basename(docking_file)
    name = file_name[:file_name.rfind(".")]
    root_dir, subdir_name = os.path.split(os.path.dirname(docking_file))
    sdf_path = os.path.join(root_dir.replace("docking", "sdf"), subdir_name, name + ".sdf")
    pdbqt_path = os.path.join(root_dir.replace("docking", "pdbqt"), subdir_name, name + ".pdbqt")

    if not os.path.exists(sdf_path):
        logger.error("Missing SDF file: %s", sdf_path)
        return []
    if not os.path.exists(pdbqt_path):
        logger.error("Missing PDBQT file: %s", pdbqt_path)
        return []

    fragments = []
    try:
        docking_result = docking_evaluate.load_docking_result(
            pdbqt_path, sdf_path, docking_file)
    except IndexError:
        logger.exception("Can't load result for: %s", docking_file)
        return []

    for model_index, molecule in enumerate(docking_result):
        discovered_fragment = docking_evaluate.extract_fragments(
            fragment_definition, molecule)
        for frag_index, atoms in enumerate(discovered_fragment):
            fragment = _create_fragment(
                fragment_definition, model_index, frag_index, molecule, atoms,
                mapping_ref, name, subdir_name)
            fragments.append(docking_evaluate.FragmentContainer(
                fragment, molecule, atoms, fragment_definition["names"]))
    return fragments


def _create_fragment(
        fragment_definition, model_index, frag_index,
        molecule, atoms, mapping_ref, name, subdir_name):
    fragment = Fragment()
    fragment.mapping = mapping_ref
    fragment.substrate = None
    fragment.name = name
    if subdir_name is not None:
        fragment.subdir_name = subdir_name
    fragment.model_index = model_index
    fragment.fragment_index = frag_index
    fragment.vina_score = molecule.prop("vina")
    fragment.atoms = []
    for atom, atom_name in zip(atoms, fragment_definition["names"]):
        fragment_atom = data_model.Fragment.Atom()
        fragment_atom.x = atom.x
        fragment_atom.y = atom.y
        fragment_atom.z = atom.z
        fragment_atom.serial = atom.serial
        fragment_atom.name = atom_name
        fragment.atoms.append(fragment_atom)
    return fragment


@functools.lru_cache(maxsize=1)
def load_protein(database, mapping_ref, conformer):
    mapping = database.read(MappingEntry, mapping_ref)
    mcsa = database.read(McsaEntry, mapping.mcsa().first().value())

    protein_pdbqt = file_formats.load_pdbqt_file(paths.target_conf_protein(conformer))

    active_site = active_site_service.active_site_from_mcsa(mcsa)
    atoms_by_names = active_site_service.reacting_atoms(
        protein_pdbqt, mapping, active_site)

    reactive_atoms = [atom for atoms in atoms_by_names.values()
                      for atom in atoms]

    return {
        "reactive_by_name": atoms_by_names,
        "reactive": reactive_atoms
    }


def _merge_dockings(files):
    new_files = set()
    dockpaths = set()
    for dfile in files:
        dockpaths.add(os.path.dirname(dfile))
    dockpaths = list(dockpaths)
    assert len(dockpaths) == 1
    files = [os.path.basename(f) for f in files if f.endswith("pdbqt")]
    for dfile in files:
        for letter in TASK["candidates"]["letters"]:
            dfile = dfile.replace(letter, "")
        new_files.add(dfile)
    new_files = [os.path.join(dockpaths[0], f) for f in new_files]
    for docking_file in new_files:
        docked_stream = []
        models = 0
        for letter in TASK["candidates"]["letters"]:
            docking_file_letter = docking_file[:-6] + letter + ".pdbqt"
            lines, models = utils.pdbqt_to_list(docking_file_letter, models)
            docked_stream.extend(lines)
        utils.list_to_pdbqt(docked_stream, docking_file)
    return new_files


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
