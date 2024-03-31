#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import json
import logging
import time
import multiprocessing as mp

from libs.configuration import initialize_configuration, TASK, ENV
from libs.remote_docking_prepare import prepare_docking_batches
from libs import utils
from libs import paths
from libs import docking
from libs import docking_evaluate
from libs import data_model
from libs.query_pubchem import query_pubchem_for_substrates
from libs.data_model import MappingEntry, McsaEntry
from libs import active_site as active_site_lib
from external import openbabel as \
    openbabel_service, mgltools as mgltools_service

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration()
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        prepare_candidates()
        logger.info("Execution time: %d s", time.time() - start_time)


def prepare_candidates():
    task_info = utils.read_task_info()
    database = data_model.create_read_only_from_directory(
        paths.knowledge_base_path())
    mcsa = database.read(McsaEntry, task_info["mcsa"])
    logger.info("Preparing candidate substrates ...")
    conformers = 1
    _export_candidates_info(database, mcsa, conformers)

    for c in range(conformers):
        if utils.is_candidate_processing_status("initial", c):
            if os.path.exists(TASK["candidates"]["prepare"]["manual"]):
                logger.info("Preparing candidates manually ...")
                _download_candidates_manually(c)
                utils.update_candidate_processing_status("sdf-downloaded", c)
            else:
                # TODO: not working yet
                _query_pubchem_for_substrates()
                utils.update_candidate_processing_status("sdf-downloaded")

        if utils.is_candidate_processing_status("sdf-downloaded", c):
            _prepare_candidates_pdbqts(c)
            utils.update_candidate_processing_status("pdbqt-prepared", c)

        if not utils.among_candidate_processing_statuses(["initial", "sdf-downloaded"], c):
            if TASK["candidates"]["docking_mode"] == "remote":
                prepare_docking_batches()
            elif TASK["candidates"]["docking_mode"] == "manual":
                # TODO
                pass

        utils.update_candidate_processing_status("docking-prepared", c)


def _export_candidates_info(database, mcsa, conformers):
    for c in range(conformers):
        pdb_id = "{}_c{}".format(TASK["knowledge_base"]["pdb"], c)
        mapping = _select_candidates_mapping(database, mcsa, pdb_id)
        name = paths.file_name_for_mapping(mapping, mcsa)
        structure_path = os.path.join(paths.proteins(), name + ".pdbqt")
        shutil.copy(structure_path, os.path.join(paths.candidates(), "protein_c{}.pdbqt".format(c)))
        _create_docking_template(mapping, mcsa, c)
    
    for c in range(conformers):
        pdb_id = "{}_c{}".format(TASK["knowledge_base"]["pdb"], c)
        mapping = _select_candidates_mapping(database, mcsa, pdb_id)
        path = os.path.join(paths.candidates(), "info_c{}.json".format(c))
        if not os.path.exists(path):
            with open(path, "w") as out_stream:
                json.dump({
                    "pdb": TASK["knowledge_base"]["pdb"],
                    "mapping": mapping.ref(),
                    "cost_function": docking_evaluate.compute_docking_cost(),
                    "status": "initial"
                }, out_stream, indent=2)


def _download_candidates_manually(conformer):
    # from libs.database import pubchempy
    import shutil
    from libs.database import pubchem
    with open(TASK["candidates"]["prepare"]["manual"], "r") as infile:
        candidate_cids = []
        for line in infile:
            try:
                int(line)
                candidate_cids.append(line.strip())
            except ValueError as ve:
                logger.error("Candidate cid:", line.strip(), "invalid")
                print(ve)
    sdf_path = paths.candidates_conf_sdf(conformer)
    os.makedirs(os.path.join(sdf_path, "subdir0001"), exist_ok=True)
    for cid in candidate_cids:
        if conformer == 0:
            cid_path = pubchem.download_3d_sdf_as_path(cid)
            shutil.copy(cid_path, os.path.join(sdf_path, "subdir0001", cid+".sdf"))
        else:
            cid_path = os.path.join(paths.candidates_conf_sdf(0), "subdir0001", cid+".sdf")
            os.symlink(cid_path, os.path.join(sdf_path, "subdir0001", cid+".sdf"))


def _query_pubchem_for_substrates():
    task_info = utils.read_task_info()
    query_pubchem_for_substrates({
        "canSMILES_value": {
            "type": "structure",
            "smiles": task_info["smiles"],
            "fragment": task_info["fragment"],
            "substrates": paths.substrates()
        },
        "XLogP3_value": {
            "type": "range",
            "value": TASK["candidates"]["prepare"]["logP"]
        },
        "MolecularWeight_value": {
            "type": "range",
            "value": TASK["candidates"]["prepare"]["MW"]
        },
        "RotatableBond_value": {
            "type": "range",
            "value": TASK["candidates"]["prepare"]["TORSDOF"]
        }
    })
    logger.info("Running with fragment ... done")


def _select_candidates_mapping(database, mcsa, pdb_id):
    mcsa_ref = mcsa.ref()
    mappings = []
    for mapping in database.iter(MappingEntry):
        if not mapping.mcsa().first().value() == mcsa_ref:
            continue
        if not mapping.source_pdb().first().value() == pdb_id:
            continue
        mappings.append(mapping)
    if len(mappings) == 0:
        raise RuntimeError("Can not find mapping for PDB: " + pdb_id +
                           " and MCSA: " + mcsa.mcsa_id())

    if TASK["knowledge_base"]["protein_chain"] is None:
        mappings.sort(key=lambda m: m.source_chain().first().value())
        return mappings[0]

    required_chain = TASK["knowledge_base"]["protein_chain"].upper()
    for mapping in mappings:
        chain = mapping.source_chain().first().value()
        if required_chain == chain:
            return mapping

    chains = [m.source_chain().first().value() for m in mappings]
    raise RuntimeError("Can not find mapping with required chain: " +
                       required_chain + "\nAvailable chains: " +
                       ", ".join(chains))


def _create_docking_template(mapping, mcsa, conformer):
    center, box = active_site_lib.compute_docking_box(mapping, mcsa)
    num_cpu_to_use = TASK["docking"]["cpu"]

    if TASK["candidates"]["docking_mode"] == "remote":
        receptor = os.path.join("in", "protein.pdbqt")
        vina_file = paths.candidates_vina_template()
    elif TASK["candidates"]["docking_mode"] == "manual":
        receptor = os.path.join(paths.candidates(), "protein_c{}.pdbqt".format(conformer))
        vina_file = os.path.join(paths.candidates(), "vina_template_c{}.txt".format(conformer))
    else:
        receptor = paths.target_protein()
        vina_file = paths.candidates_vina_template()

    docking.create_docking_file(receptor, vina_file, center, box, num_cpu_to_use)


def _prepare_candidates_pdbqts(conformer):
    if conformer == 0:
        pool = mp.Pool(processes=ENV["num_parallel_cpu"])
        logger.info("Processing candidates in %d parallel processes", ENV["num_parallel_cpu"])

        processing = []

        for files in utils.split_files_in_subdirs_to_sets(paths.candidates_conf_sdf(conformer)):
            processing.append(pool.apply_async(_prepare_candidates_set, args=(files, conformer)))

        for p in processing:
            logger.info(p.get())
        pool.close()

    else:
        prev_path = os.path.join(paths.candidates_conf_pdbqt(0), "subdir0001")
        prev_pdbqts = [os.path.join(prev_path, f) for f in os.listdir(prev_path)]
        os.makedirs(os.path.join(paths.candidates_conf_pdbqt(conformer), "subdir0001"), exist_ok=True)
        for pdbqt in prev_pdbqts:
            _pdbqt = pdbqt.replace("conformer_c0", "conformer_c{}".format(conformer))
            os.symlink(pdbqt, _pdbqt)


def _prepare_candidates_set(files, conformer):
    counter = 0
    for path_sdf in files:
        subdir_name, file_name = paths.get_subdir_filename_from_filepath(path_sdf)
        path_mol2 = os.path.join(paths.temp(), file_name.replace(".sdf", ".mol2"))
        out_dir = os.path.join(paths.candidates_conf_pdbqt(conformer), subdir_name)
        os.makedirs(out_dir, exist_ok=True)
        path_pdbqt = os.path.join(out_dir, file_name.replace(".sdf", ".pdbqt"))

        if os.path.exists(path_pdbqt):
            counter += 1
            continue
        logger.debug("Preparing candidate: %s", file_name)
        openbabel_service.convert(path_sdf, path_mol2, "sdf", "mol2")
        try:
            mgltools_service.prepare_ligand(path_mol2, path_pdbqt)
        except RuntimeError:
            logger.exception("Failed to prepare pdbqt for '{}'".format(file_name))
            continue
        finally:
            utils.remove_existing(path_mol2)
        counter += 1

    return "Processed {} files.".format(counter)


if __name__ == "__main__":
    main()
