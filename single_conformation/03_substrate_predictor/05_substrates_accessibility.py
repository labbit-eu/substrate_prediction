#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import time
import functools

from libs import file_formats
from libs import paths
from libs.data_model import McsaEntry, MappingEntry, SubstrateEntry
from libs import active_site as active_site_service
from libs import utils
from libs.configuration import initialize_configuration, TASK, ENV
from libs import accessibility1 as accessibility

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        substrates_accessibility()
        logger.info("Execution time: %d s", time.time() - start_time)


def substrates_accessibility():
    utils.transform_files(
        paths.substrates_clusters(),
        accessibility.add_accessibility_to_file,
        TASK["substrates"]["accessibility"]["input"],
        TASK["substrates"]["accessibility"]["backup"])


# region accessibility module updates

def prepare_sdf_for_protonation(database, fragments, directory):
    os.makedirs(directory, exist_ok=True)
    for fragment in fragments:
        substrate = database.read(SubstrateEntry, fragment["ref"]["substrate"])
        name = "{}-{}".format(fragment["name"], fragment["model_index"])
        output_path = os.path.join(directory, name + "_h.sdf")
        if os.path.exists(output_path):
            continue
        _prepare_sdf_file_for_protonation(fragment, substrate, name, directory)


accessibility.prepare_sdf_for_protonation = prepare_sdf_for_protonation


def _prepare_sdf_file_for_protonation(fragment, substrate, name, hydrogen_dir):
    models = accessibility.load_molecule(fragment["name"])
    molecule = models[fragment["model_index"]]

    sdf_path = paths.substrate_sdf_path(substrate)
    sdf = file_formats.load_sdf_file(sdf_path)

    pdbqt_path = paths.substrate_pdbqt_path(substrate)
    sdf_as_pdbqt = file_formats.load_pdbqt_file(pdbqt_path)

    # Update positions in SDF file to match results of docking.
    file_formats.update_sdf_positions_from_pdbqt(sdf, sdf_as_pdbqt, molecule)
    # And remove hydrogens.
    file_formats.remove_atoms_by_types(sdf, ["H"])

    temp_sdf_path = os.path.join(hydrogen_dir, name + ".sdf")
    file_formats.save_to_file(temp_sdf_path, sdf)


@functools.lru_cache(maxsize=16)
def _load_protein_atoms(database, mapping_ref):
    mapping = database.read(MappingEntry, mapping_ref)
    mcsa = database.read(McsaEntry, mapping.mcsa().first().value())

    protein_path = paths.protein_pdbqt_path(mapping, mcsa)
    protein_pdbqt = file_formats.load_pdbqt_file(protein_path)

    active_site = active_site_service.active_site_from_mcsa(mcsa)
    atoms_by_names = active_site_service.reacting_atoms(
        protein_pdbqt, mapping, active_site)

    return {
        "reactive_by_name": atoms_by_names,
        "all": protein_pdbqt.atoms(),
    }


accessibility.load_protein_atoms = _load_protein_atoms


@functools.lru_cache(maxsize=64)
def _load_molecule(name):
    path = paths.complexes_pdbqt_path(name)
    return file_formats.load_pdbqt_models_file(path)


accessibility.load_molecule = _load_molecule

# endregion

if __name__ == "__main__":
    main()
