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


LEPROT = "PATH_TO_LEPRO_LEDOCK"
LEDOCK = "PATH_TO_LEDOCK_BINARY"
OBABEL = "PATH_TO_OPENBABEL"

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

def count_conformers(pdb_path, out_path):
    pdbs = [p for p in os.listdir(pdb_path) if p.endswith("pdb")]
    return len(pdbs)

LEDOCK_CONFIG = """Receptor
{RECEPTOR}

RMSD
0.5

Binding pocket
{XMIN} {XMAX}
{YMIN} {YMAX}
{ZMIN} {ZMAX}

Number of binding poses
200

Ligands list
{LIGSPATH}

END
"""

def _prepare_substrate_mol2(sdfin, mol2out):
    if os.path.isfile(mol2out):
        return
    try:
        command = [OBABEL, "-i", "sdf", sdfin, "-o", "mol2", "-O", mol2out]
        _execute_command(command)
    except RuntimeError:
        print("Failed to prepare mol2 for '{}'".format(sdfin))

def _prepare_receptor(pdbin, pdbout):
    if os.path.isfile(pdbout):
        return
    if os.path.isfile("prot.pdb"):
        raise RuntimeError("Temporary PDB file already present, remove it first.")
    try:
        command = [LEPROT, pdbin]
        _execute_command(command)
        os.rename("pro.pdb", pdbout)
        os.remove("dock.in")
        assert not os.path.isfile("pro.pdb")
        assert not os.path.isfile("dock.in")
    except RuntimeError:
        print("Failed to prepare receptor file '{}'".format(pdbin))

def prepare_docking_files(pdbin, ligslist, ligspath, box, outpath):
    center = box["center"]
    size = box["size"]
    xmin, xmax = center[0]-(size[0]/2.0), center[0]+(size[0]/2.0)
    ymin, ymax = center[1]-(size[1]/2.0), center[1]+(size[1]/2.0)
    zmin, zmax = center[2]-(size[2]/2.0), center[2]+(size[2]/2.0)
    config = LEDOCK_CONFIG.format(RECEPTOR=pdbin, XMIN=xmin, XMAX=xmax,
                                  YMIN=ymin, YMAX=ymax, ZMIN=zmin, ZMAX=zmax,
                                  LIGSPATH=ligspath)
    with open(ligspath, "w") as fs:
        fs.write(ligslist)
    with open(outpath, "w") as outstream:
        outstream.write(config)

def prepare_ledocks(database):
    box = database["DOCKING_BOX"]
    os.makedirs("receptors", exist_ok=True)
    for uniprot in database["UNIPROT"].values():
        # Prepare mol2 files for substrates and make the ligands list
        ligslist = ""
        for substrate in uniprot["substrates"]:
            if not os.path.isfile("substrates/{}.mol2".format(substrate)):
                _prepare_substrate_mol2("substrates/{}.sdf".format(substrate),
                                        "substrates/{}.mol2".format(substrate))
            ligslist += "substrates/{}.mol2\n".format(substrate)
        # Prepare the pdb files for receptors and make the config files for docking
        for pdbid in uniprot["pdbs"]:
            os.makedirs(pdbid, exist_ok=True)
            os.makedirs("{}/ledock".format(pdbid), exist_ok=True)
            conformers = count_conformers("../01_fpocket/{}".format(pdbid), "receptors")
            for i in range(conformers):
                pdb = "receptors/{}_c{}.pdb".format(pdbid, i)
                _pdb = "receptors/{}_c{}_LD.pdb".format(pdbid, i)
                if not os.path.isfile(_pdb):
                    _prepare_receptor(pdb, _pdb)
                prepare_docking_files(_pdb, ligslist, "{}/ledock/ligands.txt".format(pdbid),
                                      box, "{}/ledock/{}_c{}_LD.in".format(pdbid, pdbid, i))

def run_ledocks(database):
    for uniprot in database["UNIPROT"].values():
        for pdbid in uniprot["pdbs"]:
            conformers = count_conformers("../01_fpocket/{}".format(pdbid), "receptors")
            for i in range(conformers):
                config_file = "{}/ledock/{}_c{}_LD.in".format(pdbid, pdbid, i)
                try:
                    command = [LEDOCK, config_file]
                    _execute_command(command)
                except RuntimeError:
                    print("Failed to run docking for '{}_c{}'".format(pdbid, i))
                for substrate in uniprot["substrates"]:
                    assert os.path.isfile("substrates/{}.dok".format(substrate)),\
                           "Failed to find {} at {} with pdb {}".format(substrate, uniprot, pdbid)
                    os.rename("substrates/{}.dok".format(substrate),
                              "{}/ledock/{}_c{}_@_{}.dok".format(pdbid, pdbid, i, substrate))
                    assert not os.path.isfile("substrates/{}.dok".format(substrate)),\
                           "Failed to remove {} at {} with pdb {}".format(substrate, uniprot, pdbid)

def main(datafile):
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    prepare_ledocks(database)
    run_ledocks(database)


# python3 04_run_ledock.py json_database.json
main(sys.argv[1])
