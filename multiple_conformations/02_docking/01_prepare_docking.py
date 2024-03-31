#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  01_prepare_docking.py
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import os
import sys
import json
import subprocess


PYMOL_PYTHON = "PATH_TO_PYMOL_PYTHON"
MGL_PYTHON = "PATH_TO_PYTHONSH_MGLTOOLS"
REF_PDB = "REFERENCE_PDB_FILE"
PREPARE_RECEPTOR4 = "PATH_TO_MGLTOOLS_PREPARE_RECEPTOR4"

def _execute_command(command, shell=False, cwd=None):
    with subprocess.Popen(
            command,
            cwd=cwd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE) as process:
        out, err = process.communicate()
        errcode = process.returncode
    if errcode == 0:
        return
    message = "Shell command failed: \n"
    message += " ".join(command)
    message += "\nerror code:"
    message += str(errcode)
    if out is not None:
        message += "\nstdout:\n"
        message += out.decode("utf-8")
    if err is not None:
        message += "\nstderr:\n"
        message += err.decode("utf-8")
    raise RuntimeError(message)

def _pymol_align_script(ref, inpdb, outpdb):
    return """#! {}
import pymol
from pymol import cmd
pymol.pymol_argv = ["pymol", "-qc"]
pymol.finish_launching()
ref = "reference"
pdb = "protein"
cmd.load("{}", ref)
cmd.load("{}", pdb)
results = cmd.align(pdb, ref)
cmd.remove("hydrogens")
cmd.save("{}", pdb)
""".format(PYMOL_PYTHON, ref, inpdb, outpdb)

def align_receptors(pdb_path, out_path):
    pdbs = [p for p in os.listdir(pdb_path) if p.endswith("pdb")]
    for pdb in pdbs:
        aligned_pdb = "{}/{}".format(out_path, pdb)
        if not os.path.isfile(aligned_pdb):
            pymol_script = _pymol_align_script(REF_PDB, "{}/{}".format(pdb_path, pdb), aligned_pdb)
            script_path = "tmp_align.py"
            with open(script_path, "w") as temp_stream:
                temp_stream.write(pymol_script)
            command = " ".join([PYMOL_PYTHON, script_path])
            _execute_command(command, shell=True)
            os.system("sed -i 's/Oa/ O/g' receptors/{}".format(pdb))
            os.system("sed -i 's/Op/ O/g' receptors/{}".format(pdb))
            os.remove(script_path)
    return len(pdbs)

VINA_CONFIG = """center_x = {CENTER_X:f}
center_y = {CENTER_Y:f}
center_z = {CENTER_Z:f}
size_x = {SIZE_X:f}
size_y = {SIZE_Y:f}
size_z = {SIZE_Z:f}
cpu = {CPU:d}
exhaustiveness = 10
num_modes = 200
energy_range = 10
"""

def _prepare_protein_pdbqt(pdbin, pdbqtout):
    try:
        args = ["-A", "checkhydrogens", "-U", "nphs_lps"]
        command = [MGL_PYTHON, PREPARE_RECEPTOR4] + ["-r", pdbin, "-o", pdbqtout] + args
        _execute_command(command)
    except RuntimeError:
        print("Failed to prepare pdbqt for '{}'".format(pdbin))
        return None
    return pdbqtout

def prepare_docking_files(pdbid, substrate, box, conformers=6):
    center = box["center"]
    size = box["size"]
    if not os.path.isfile("substrates/{}.pdbqt".format(substrate)):
        raise RuntimeError("Substrates file does not exist", substrate)
    for i in range(conformers):
        pdb = "receptors/{}_c{}.pdb".format(pdbid, i)
        pdbqt = "receptors/{}_c{}.pdbqt".format(pdbid, i)
        if not os.path.isfile(pdbqt):
            _prepare_protein_pdbqt(pdb, pdbqt)
        config = VINA_CONFIG.format(CENTER_X=center[0], CENTER_Y=center[1],
                                    CENTER_Z=center[2], SIZE_X=size[0],
                                    SIZE_Y=size[1], SIZE_Z=size[2], CPU=10)
        config += "receptor = {}\n".format(pdbqt)
        config += "ligand = substrates/{}.pdbqt\n".format(substrate)
        config += "out = {}/outputs/{}_c{}_@_{}.pdbqt\n".format(pdbid, pdbid, i, substrate)
        config += "log = {}/logs/{}_c{}_@_{}.log\n".format(pdbid, pdbid, i, substrate)
        config_path = "{}/configs/{}_c{}_@_{}.txt".format(pdbid, pdbid, i, substrate)
        with open(config_path, "w") as outstream:
            outstream.write(config)

def main(datafile):
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    box = database["DOCKING_BOX"]
    os.makedirs("receptors", exist_ok=True)
    for uniprot in database["UNIPROT"].values():
        for pdb in uniprot["pdbs"]:
            os.makedirs(pdb, exist_ok=True)
            os.makedirs("{}/configs".format(pdb), exist_ok=True)
            os.makedirs("{}/logs".format(pdb), exist_ok=True)
            os.makedirs("{}/outputs".format(pdb), exist_ok=True)
            conformers = align_receptors("../03_fpocket/{}".format(pdb), "receptors")
            for substrate in uniprot["substrates"]:
                prepare_docking_files(pdb, substrate, box, conformers=conformers)
    conformers = align_receptors("../03_fpocket/{}".format(database["QUERY_PROTEIN"]["PDB"]), "receptors")
    for i in range(conformers):
        pdb = "receptors/{}_c{}.pdb".format(database["QUERY_PROTEIN"]["PDB"], i)
        pdbqt = "receptors/{}_c{}.pdbqt".format(database["QUERY_PROTEIN"]["PDB"], i)
        if not os.path.isfile(pdbqt):
            _prepare_protein_pdbqt(pdb, pdbqt)


# python3 01_prepare_docking.py json_database.json
main(sys.argv[1])
