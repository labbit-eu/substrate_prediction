#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import logging

from libs import utils
from libs.remote_docking import candidates_remote_docking
from libs import paths
from libs import docking
from libs.configuration import initialize_configuration, TASK, ENV, REMOTE

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
            if utils.is_candidate_processing_status("docking-completed", c):
                logger.warning("Screening of candidates has already finished for this task. Proceed with further analysis")
                exit(0)
            if not utils.among_candidate_processing_statuses(["docking-prepared", "ready_for_transfer2remote",
                                                              "transferred2remote", "docking-failed"], c):
                logger.error("Candidates were not prepared for docking. Please run EXT08_candidates_prepare.py first")
                exit(1)
            else:
                if TASK["candidates"]["docking_mode"] == "remote":
                    if utils.is_candidate_processing_status("docking-failed") and not REMOTE["update_template_file"]:
                        logger.error("Previous screening with the current setting failed, please modify the template file "
                                     "and any relevant settings and set the parameter \"update_template_file\" in [REMOTE] "
                                     "section of the configuration file to True to proceed with a new screening campaign.")
                        exit(1)
                    candidates_remote_docking()
                elif TASK["candidates"]["docking_mode"] == "manual":
                    if utils.is_candidate_processing_status("docking-failed", c):
                        logger.error("Previous screening with the current setting failed, please modify the template file "
                                     "and any relevant settings and set the parameter \"update_template_file\" in [REMOTE] "
                                     "section of the configuration file to True to proceed with a new screening campaign.")
                        exit(1)
                    candidates_docking_manual(c)
                else:
                    candidates_docking()

    logger.info("Execution time: %d s", time.time() - start_time)


def candidates_docking_manual(conformer):
    root_dir = paths.candidates_conf_pdbqt(conformer)
    num_files = 0
    for subdir_name in os.listdir(root_dir):
        num_files += len(os.listdir(os.path.join(root_dir, subdir_name)))
    num_files *= len(TASK["candidates"]["letters"])
    logger.info("Docking: %d candidate substrates in conformer %d", num_files, conformer)
    
    num_lig = 0
    for letter in TASK["candidates"]["letters"]:
        for subdir_name in os.listdir(root_dir):
            output_dir = os.path.join(paths.candidates_conf_docking(conformer), subdir_name)
            os.makedirs(output_dir, exist_ok=True)
            vina_file = os.path.join(paths.candidates(), "vina_template_c{}.txt".format(conformer))
            for file_name in os.listdir(os.path.join(root_dir, subdir_name)):
                output_file = "{}{}".format(file_name[:-6], letter)
                ligand_path = os.path.join(root_dir, subdir_name, file_name)
                num_lig += 1
                logger.info("Docking: %d/%d file", num_lig,  num_files)
                docking.dock_ligand_vina(vina_file, output_dir, output_file, ligand_path)
    
    utils.update_candidate_processing_status("docking-completed", conformer)


def candidates_docking():
    root_dir = paths.candidates_pdbqt()
    num_files = 0
    for subdir_name in os.listdir(root_dir):
        num_files += len(os.listdir(os.path.join(root_dir, subdir_name)))

    logger.info("Docking: %d candidate substrates", num_files)
    
    num_lig = 0
    for subdir_name in os.listdir(root_dir):
        output_dir = os.path.join(paths.candidates_docking(), subdir_name)
        os.makedirs(output_dir, exist_ok=True)
        for file_name in os.listdir(os.path.join(root_dir, subdir_name)):
            ligand_path = os.path.join(root_dir, subdir_name, file_name)
            num_lig += 1
            logger.info("Docking: %d/%d file", num_lig,  num_files)
            docking.dock_ligand_vina(paths.candidates_vina_template(), output_dir, file_name[:-6], ligand_path)

    utils.update_candidate_processing_status("docking-completed")


if __name__ == "__main__":
    main()

