#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import time
import os
import functools

from libs import file_formats
from libs.data_model import McsaEntry, MappingEntry
from libs import active_site as active_site_service
from libs import paths
from libs import utils
from libs import accessibility2 as accessibility
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        candidates_accessibility()
        logger.info("Execution time: %d s", time.time() - start_time)


def candidates_accessibility():
    utils.transform_files(
        paths.candidates_filtered(),
        accessibility.add_accessibility_to_file,
        TASK["candidates"]["accessibility"]["input"],
        TASK["candidates"]["accessibility"]["backup"])


# region accessibility module updates
accessibility.MODE = "candidates"


@functools.lru_cache(maxsize=16)
def _load_protein_atoms(database, mapping_ref, conformer):
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


accessibility.load_protein_atoms = _load_protein_atoms


@functools.lru_cache(maxsize=64)
def _load_molecule(fragment_name, subdir_name, conformer):
    path = get_fragment_path("docking", fragment_name, subdir_name, conformer)
    return file_formats.load_pdbqt_models_file(path)


accessibility.load_molecule = _load_molecule


def get_fragment_path(type_, fragment_name, subdir_name, conformer):
    root_dir = paths.candidates()
    if type_ == "docking":
        return os.path.join(root_dir, "conformer_c{}".format(conformer), "docking", subdir_name, fragment_name + ".pdbqt")
    elif type_ == "pdbqt":
        return os.path.join(root_dir, "conformer_c{}".format(conformer), "pdbqt", subdir_name, fragment_name + ".pdbqt")
    elif type_ == "sdf":
        return os.path.join(root_dir, "conformer_c{}".format(conformer), "sdf", subdir_name, fragment_name + ".sdf")
    else:
        raise RuntimeError("Invalid path type " + type_)


# def prepare_sdf_for_protonation(database, fragments, directory):
def prepare_sdf_for_protonation(database, fragments):
    for fragment in fragments:
        _conf = fragment["ref"]["mapping"].split(",")[1]
        _conf = _conf.split("-")[0]
        _conf = _conf.split("_")[-1][1:]
        conformer = int(_conf)
        directory = os.path.join(paths.temp(), accessibility._get_temp_directory(conformer))
        os.makedirs(directory, exist_ok=True)
        name = "{}-{}".format(fragment["name"], fragment["model_index"])
        os.makedirs(os.path.join(directory, fragment["subdir_name"]), exist_ok=True)
        output_path = os.path.join(directory, fragment["subdir_name"], name + "_h.sdf")
        if os.path.exists(output_path):
            continue
        _prepare_sdf_file_for_protonation(fragment, name, directory, conformer)


accessibility.prepare_sdf_for_protonation = prepare_sdf_for_protonation


def _prepare_sdf_file_for_protonation(fragment, name, hydrogen_dir, conformer):
    models = accessibility.load_molecule(fragment["name"], fragment["subdir_name"], conformer)
    molecule = models[fragment["model_index"]]

    sdf_path = get_fragment_path("sdf", fragment["name"], fragment["subdir_name"], conformer)
    sdf = file_formats.load_sdf_file(sdf_path)

    pdbqt_path = get_fragment_path("pdbqt", fragment["name"], fragment["subdir_name"], conformer)
    sdf_as_pdbqt = file_formats.load_pdbqt_file(pdbqt_path)

    # Update positions in SDF file to match results of docking.
    file_formats.update_sdf_positions_from_pdbqt(sdf, sdf_as_pdbqt, molecule)
    # And remove hydrogens.
    file_formats.remove_atoms_by_types(sdf, ["H"])

    temp_sdf_path = os.path.join(hydrogen_dir,  fragment["subdir_name"], name + ".sdf")
    file_formats.save_to_file(temp_sdf_path, sdf)


# endregion

if __name__ == "__main__":
    main()
