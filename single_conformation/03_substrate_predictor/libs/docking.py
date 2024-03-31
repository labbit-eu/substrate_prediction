#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging

from external import autodock_vina
from libs.configuration import TASK

logger = logging.getLogger(__name__)


def get_torsdof_from_pdbqt(pdbqt_path):
    with open(pdbqt_path) as pdbqt_stream:
        for line in pdbqt_stream.readlines():
            if "TORSDOF" in line:
                return line.split()[1]
    logger.error("TORSDOF not defined for file: %s", pdbqt_path)
    return None


def dock_ligand_vina(configuration_path, output_dir, ligand_filename_without_ext, input_pdbqt_file=None):
    """
    presumes that docked pdbqt file and log files are in output_dir and both differ by their file extension only
    => they share ligand_filename_without_ext
    """
    output_path_pdbqt = os.path.join(output_dir, ligand_filename_without_ext + ".pdbqt")
    output_path_log = os.path.join(output_dir, ligand_filename_without_ext + ".log")
        
    if test_vina_results(output_path_log, output_path_pdbqt):
        return
    try:
        autodock_vina.execute(configuration_path, [output_path_pdbqt, output_path_log], input_pdbqt_file)
    except RuntimeError:
        logger.exception("Failed to dock file {}".format(input_pdbqt_file))


def test_vina_results(output_path_log, output_path_pdbqt):
    """
    test completeness and consistency of files with docking results
    """
    if os.path.exists(output_path_log) and os.path.exists(output_path_pdbqt):
        with open(output_path_log, "r") as log_stream, open(output_path_pdbqt, "r") as pdbqt_stream:
            log_lines = log_stream.readlines()
            pdbqt_lines = pdbqt_stream.readlines()
            if ("Writing output ... done." in log_lines[-1]) and ("ENDMDL" in pdbqt_lines[-1]):
                num_modes_log = int(log_lines[-2].split()[0])
                num_modes_pdbqt = max([int(x.split()[1]) for x in pdbqt_lines if "MODEL" in x])
                if num_modes_log == num_modes_pdbqt:
                    return True
    return False


def create_docking_file(receptor_path, config_out_path, center, box, num_cpu, out_name=None, ligand_path=None):
    if num_cpu is None:
        num_cpu = TASK["docking"]["cpu"]

    with open(config_out_path, "w") as output_stream:
        output_stream.write("""center_x = {}
center_y = {}
center_z = {}

size_x = {}
size_y = {}
size_z = {}

cpu = {}
exhaustiveness = {}
#seed = {}
num_modes = {}
energy_range = {}

receptor = {}
""".format(
            center[0], center[1], center[2],
            box[0], box[1], box[2],
            num_cpu,
            TASK["docking"]["exhaustiveness"],
            TASK["docking"]["seed"],
            TASK["docking"]["num_modes"],
            TASK["docking"]["energy_range"],
            receptor_path
        ))
        if ligand_path is not None:
            output_stream.write("ligand = {}\n".format(ligand_path))

        if out_name is not None:
            output_stream.write("""out = {}
log = {}
""".format(out_name + ".pdbqt", out_name + ".log",))

