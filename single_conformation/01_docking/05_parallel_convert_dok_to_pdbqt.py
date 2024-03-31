#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  05_parallel_convert_dok_to_pdbqt.py
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carse@amu.edu.pl>
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
import tempfile
import subprocess
import numpy as np
from typing import List, Dict
from multiprocessing import Pool
from dataclasses import dataclass, field

MGL_PYTHON = "PATH_TO_PYTHONSH_MGLTOOLS"
PREPARE_LIGAND = "PATH_TO_PREPARE_LIGAND"
np.set_printoptions(precision=3)

@dataclass(order=True)
class Atom:
    atid: int = field(compare=True, hash=True)
    name: str = field(compare=True, hash=True)
    pos: np.array = field()

@dataclass(order=True)
class PDB:
    pose: int = field(compare=True, hash=True)
    score: float = field(compare=True, hash=True)
    atoms: Dict = field(repr=True, default_factory=dict, hash=False)

@dataclass(order=True)
class PDBQT:
    pose: int = field(compare=True, hash=True)
    score: float = field(compare=True, hash=True)
    atoms: Dict = field(repr=True, default_factory=dict, hash=False)
    content: List = field(repr=False, default_factory=list, hash=False)


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

def _check_elements(content):
    for i in range(len(content)):
        if content[i][12:15] == " CL":
            content[i] = "{}CL {}".format(content[i][:12], content[i][15:])
        if content[i][12:15] == " BR":
            content[i] = "{}BR {}".format(content[i][:12], content[i][15:])

def make_first_model_pdbqt(infile, temp_path):
    with open(infile, "r") as instream:
        line = instream.readline()
        assert line.startswith("REMARK Docking time")
        output = []
        line = instream.readline()
        while not line.startswith("END"):
            output.append(line)
            line = instream.readline()
        output.append(line)
        _check_elements(output)
    tmp_pdb = os.path.join(temp_path, "tmp_1.pdb")
    tmp_pdbqt = os.path.join(temp_path, "tmp_1.pdbqt")
    with open(tmp_pdb, "w") as temp:
        temp.write("".join(output))
    command = [MGL_PYTHON, PREPARE_LIGAND, "-l", tmp_pdb, "-o", tmp_pdbqt, "-U", "\"''\""]
    try:
        _execute_command(command)
    except RuntimeError:
        print("Failed to prepare ligand with MglTools.", tmp_pdb, tmp_pdbqt)
    return tmp_pdb, tmp_pdbqt

def get_docked_poses(infile):
    docked_poses = []
    with open(infile, "r") as instream:
        i = 1
        line = instream.readline()
        assert line.startswith("REMARK Docking time")
        while line:
            if line.startswith("REMARK Cluster"):
                chunks = line.split()
                docked = PDB(i, float(chunks[-2]))
                while not line.startswith("END"):
                    if line.startswith("ATOM"):
                        chunks = line.split()
                        atom = Atom(int(chunks[1]), chunks[2], np.array([float(chunks[5]),float(chunks[6]),float(chunks[7])]))
                        docked.atoms[int(chunks[1])] = atom
                    line = instream.readline()
                docked_poses.append(docked)
                i += 1
            line = instream.readline()
    docked_poses.sort()
    return docked_poses

def map_pdbqt_pdb(pdbfile, pdbqtfile):
    with open(pdbfile, "r") as f:
        chunks = f.readline().split()
        pdb = PDB(1, float(chunks[-2]))
        for line in f:
            if line.startswith("ATOM"):
                chunks = line.split()
                atom = Atom(int(chunks[1]), chunks[2].strip(), np.array([float(chunks[5]),float(chunks[6]),float(chunks[7])]))
                pdb.atoms[int(chunks[1])] = atom
    with open(pdbqtfile, "r") as f:
        pdbqt = PDBQT(1, pdb.score)
        for line in f:
            if line.startswith("ATOM "):
                chunks = line.split()
                atom = Atom(int(chunks[1]), chunks[2].strip(), np.array([float(chunks[5]),float(chunks[6]),float(chunks[7])]))
                pdbqt.atoms[int(chunks[1])] = atom
            pdbqt.content.append(line)
    pdbqtid_to_pdbid = {}
    for qt_atom in pdbqt.atoms.values():
        for pb_atom in pdb.atoms.values():
            if sum(abs(qt_atom.pos-pb_atom.pos)) < 0.001:
                pdbqtid_to_pdbid[qt_atom.atid] = pb_atom.atid
                break
        else:
            print("No equivalent atom found for", qt_atom)
    return pdbqtid_to_pdbid, pdbqt

def map_pdbqt_mol2(pdbqtfile, mol2file):
    with open(pdbqtfile, "r") as f:
        pdbqt = PDBQT(1, 0.0)
        for line in f:
            if line.startswith("ATOM "):
                chunks = line.split()
                x = np.round(float(chunks[6]), 3)
                y = np.round(float(chunks[7]), 3)
                z = np.round(float(chunks[8]), 3)
                atom = Atom(int(chunks[1]), chunks[2].strip(), np.array([x, y, z]))
                pdbqt.atoms[int(chunks[1])] = atom
            pdbqt.content.append(line)
    with open(mol2file, "r") as f:
        mol2 = PDB(1, 0.0)
        line = f.readline()
        while line:
            if line.startswith("@<TRIPOS>ATOM"):
                line = f.readline()
                while not line.startswith("@<TRIPOS>"):
                    if len(line.strip()) > 5:
                        chunks = line.split()
                        x = np.round(float(chunks[2]), 3)
                        y = np.round(float(chunks[3]), 3)
                        z = np.round(float(chunks[4]), 3)
                        atom = Atom(int(chunks[0]), chunks[1].strip(), np.array([x, y, z]))
                        mol2.atoms[int(chunks[0])] = atom
                    line = f.readline()
                break
            line = f.readline()
    
    pdbqtid_to_mol2id = {}
    for qt_atom in pdbqt.atoms.values():
        for m2_atom in mol2.atoms.values():
            if (abs(qt_atom.pos-m2_atom.pos) < 0.01).all():
                pdbqtid_to_mol2id[qt_atom.atid] = m2_atom.atid
                break
        else:
            print("No equivalent atom found for", qt_atom, pdbqtfile, mol2file)
    return pdbqtid_to_mol2id, pdbqt

def convert_dok_to_pdbqt(poses, atom_map, pdbqt_patt, outpath):
    output = []
    for pose in poses:
        output.append("MODEL {}\n".format(pose.pose))
        output.append("REMARK VINA RESULT:{:>10.1f}{:>11.3f}{:>11.3f}\n".format(pose.score, 0, 0))
        for line in pdbqt_patt.content:
            if line.startswith("ATOM "):
                chunks = line.split()
                x, y, z = pose.atoms[atom_map[int(chunks[1])]].pos
                newline = "{}{:>8.3f}{:>8.3f}{:>8.3f}{}".format(line[:30], x, y, z, line[54:])
                output.append(newline)
            else:
                output.append(line)
        output.append("ENDMDL\n")
    with open(outpath, "w") as fout:
        fout.write("".join(output))

def convert_parallel(substrates, pdbid, conformer):
    for substrate in substrates:
        dokfile = "{}/ledock/{}_c{}_@_{}.dok".format(pdbid, pdbid, conformer, substrate)
        poses = get_docked_poses(dokfile)
        mol2file = os.path.join("substrates", "{}.mol2".format(substrate))
        pdbqtfile = os.path.join("substrates", "{}.pdbqt".format(substrate))
        atom_map, pdbqt_patt = map_pdbqt_mol2(pdbqtfile, mol2file)
        out_pdbqtfile = "{}/ledock/{}_c{}_@_{}_LD.pdbqt".format(pdbid, pdbid, conformer, substrate)
        convert_dok_to_pdbqt(poses, atom_map, pdbqt_patt, out_pdbqtfile)

def main(datafile):
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    for uniprot in database["UNIPROT"].values():
        for pdbid in uniprot["pdbs"]:
            conformers = 1
            pool = Pool(conformers)
            for i in range(conformers):
                pool.apply_async(convert_parallel, args=(uniprot["substrates"], pdbid, i))
            pool.close()
            pool.join()


if __name__ == "__main__":
    # python3 05_parallel_convert_dok_to_pdbqt.py json_database.json
    main(sys.argv[1])

