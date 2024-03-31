#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from external import process
from libs import paths
from libs.utils import remove_existing
from libs.configuration import ENV


def add_hydrogens_batch(directory, mode):
    if mode == "substrates":
        _add_hydrogens_batch(directory)
    else:
        pool = mp.Pool(processes=ENV["num_parallel_cpu"])
        processing = []
    
        for subdir_name in sorted(os.listdir(directory)):
            if subdir_name.startswith("subdir"):
                dir_to_process = os.path.join(directory, subdir_name)
                processing.append(pool.apply_async(_add_hydrogens_batch, args=(dir_to_process, subdir_name)))
        
        for p in processing:
            p.get()
        pool.close()


def _add_hydrogens_batch(directory, subdir_name=None):
    if subdir_name is not None:
        script_path = os.path.join(paths.temp(), subdir_name+"add_hydrogens_to_ligand_batch.py")
    else:
        script_path = os.path.join(paths.temp(), "add_hydrogens_to_ligand_batch.py")
    with open(script_path, "w") as output_stream:
        output_stream.write(_create_batch_script(directory))

    process.execute([ENV["python_pymol"], script_path],
                    shell=ENV["python_pymol_shell"])
    remove_existing(script_path)


def _create_batch_script(directory):
    return """
import os
    
import pymol
from pymol import cmd
pymol.pymol_argv = ["pymol", "-qc"]
pymol.finish_launching()

directory = "{}"
for file_name in os.listdir(directory):
    if file_name.endswith("_h.sdf"):
        continue
    in_path = os.path.join(directory, file_name)
    out_path = os.path.join(directory, file_name[0:-4] + "_h.sdf")
    
    cmd.load(in_path, "mol")
    cmd.h_add("mol")
    cmd.save(out_path, "mol")
    cmd.delete("*")
    os.remove(in_path)
    
    """.format(directory.replace("\\", "/"))
