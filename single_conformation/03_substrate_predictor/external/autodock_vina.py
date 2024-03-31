#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from external import process
from libs.configuration import ENV


def execute(configuration_path, outputs=None, in_ligand_path=None,):
    # outputs = [output_path_pdbqt, output_path_log]
    cwd = os.path.dirname(configuration_path)
    file_name = os.path.basename(configuration_path)
    vina_path = ENV["vina"]
    if vina_path.startswith("."):
        # Relative path from root we need to update it, to be relative
        # from cwd. which is in this case data/{}/{}
        vina_path = os.path.join("..", "..", "..", vina_path)
    command = [vina_path, "--config", file_name]
    
    if outputs is not None:
        if len(outputs) == 2: 
            command = command + ["--out", outputs[0], "--log", outputs[1]]
        elif len(outputs) == 1: 
            command = command + ["--out", outputs[0]]
        else: 
            raise RuntimeError("Error in definition of outputs parameters in autodock_vina.execute")

    if in_ligand_path is not None:
        command = command + ["--ligand", in_ligand_path]

    process.execute(command, cwd=cwd, shell=ENV["vina_shell"])
