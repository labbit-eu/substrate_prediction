#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import sys
import json
from typing import List, Dict
from dataclasses import dataclass, field


@dataclass(order=True)
class DOK:
    pose: int = field(compare=True, hash=True)
    score: float = 0.0
    content: List = field(repr=False, default_factory=list, hash=False)
    def _check_elements(self):
        for i in range(len(self.content)):
            if self.content[i][12:15] == " CL":
                self.content[i] = "{}CL {}".format(self.content[i][:12], self.content[i][15:])
            if self.content[i][12:15] == " BR":
                self.content[i] = "{}BR {}".format(self.content[i][:12], self.content[i][15:])
    def __str__(self):
        return "".join(self.content)

@dataclass(order=True)
class PDBQT:
    pose: int = field(compare=True, hash=True)
    score: float = 0.0
    atoms: Dict = field(repr=False, default_factory=dict, hash=False)
    content: List = field(repr=False, default_factory=list, hash=False)
    def check_atom_order(self, other):
        if len(self.atoms.keys()) != len(other.atoms.keys()):
            return "PDBQT objects must have the same number of atoms. Got {} and {}".format(len(self.atoms), len(other.atoms))
        for atom_id in self.atoms.keys():
            if self.atoms[atom_id] != other.atoms[atom_id]:
                return "Atoms with id {} have different positions in PDBQT files, please check model {}.\n{}".format(atom_id, other.pose, "".join(other.content))
        return None


def _count_conformers(pdb_path):
    pdbs = [p for p in os.listdir(pdb_path) if p.endswith("pdb")]
    return len(pdbs)

def get_docked_poses(infile):
    docked_poses = []
    with open(infile, "r") as instream:
        i = 1
        line = instream.readline()
        assert line.startswith("REMARK Docking time")
        while line:
            if line.startswith("REMARK Cluster"):
                chunks = line.split()
                docked = DOK(i, float(chunks[-2]))
                while not line.startswith("END"):
                    docked.content.append(line)
                    line = instream.readline()
                docked.content.append(line)
                docked._check_elements()
                docked_poses.append(docked)
                i += 1
            line = instream.readline()
    return docked_poses

def check_pdbqt_order(infile, debug=False):
    if debug:
        output = "Order to follow (MODEL 1):\n"
    else:
        output = ""
    with open(infile, "r") as instream:
        line = instream.readline()
        chunks = line.split()
        pose = int(chunks[1])
        line = instream.readline()
        chunks = line.split()
        score = float(chunks[3])
        order_pdbqt = PDBQT(pose, score)
        while not line.startswith("ENDMDL"):
            if line.startswith("ATOM"):
                chunks = line.strip().split()
                atid = int(chunks[1])
                assert atid not in order_pdbqt.atoms, "Atom {} twice in molecule {}".format(atid, pose)
                order_pdbqt.atoms[atid] = (chunks[2], chunks[-1])
                order_pdbqt.content.append(line)
            line = instream.readline()
        line = instream.readline()
        if debug:
            output += "".join(order_pdbqt.content)
            output += "\n"
        
        while line:
            if line.startswith("MODEL"):
                chunks = line.strip().split()
                pose = int(chunks[1])
                line = instream.readline()
                chunks = line.split()
                score = float(chunks[3])
                to_test = PDBQT(pose, score)
                while not line.startswith("ENDMDL"):
                    if line.startswith("ATOM"):
                        chunks = line.strip().split()
                        atid = int(chunks[1])
                        assert atid not in to_test.atoms, "Atom {} twice in molecule {}".format(atid, pose)
                        to_test.atoms[atid] = (chunks[2], chunks[-1])
                        to_test.content.append(line)
                    line = instream.readline()
                msg = order_pdbqt.check_atom_order(to_test)
                if msg is not None:
                    output += msg
                to_test = None
            line = instream.readline()
    return output

def make_poses(infile):
    docked_poses = get_docked_poses(infile)
    docked_poses.sort()
    print(docked_poses)
    for i, dp in enumerate(docked_poses):
        with open("test_{}.pdb".format(i+1), "w") as temp:
            temp.write(str(dp))

def get_pdbqt_models(infile):
    models = []
    with open(infile, "r") as instream:
        line = instream.readline()
        while line:
            if line.startswith("MODEL"):
                chunks = line.strip().split()
                pose = int(chunks[1])
                line = instream.readline()
                chunks = line.split()
                score = float(chunks[3])
                _pdbqt = PDBQT(pose, score)
                while not line.startswith("ENDMDL"):
                    if line.startswith("ATOM"):
                        chunks = line.strip().split()
                        atid = int(chunks[1])
                        assert atid not in _pdbqt.atoms, "Atom {} twice in molecule {}".format(atid, pose)
                        _pdbqt.atoms[atid] = (chunks[2], chunks[-1])
                    _pdbqt.content.append(line)
                    line = instream.readline()
                _pdbqt.content.append(line)
                models.append(_pdbqt)
            line = instream.readline()
    return models

def check_order_intra(database):
    output = ""
    for uniprot in database["UNIPROT"].values():
        for pdbid in uniprot["pdbs"]:
            files = [f for f in os.listdir(os.path.join(pdbid, "ledock")) if f.endswith(".pdbqt")]
            for f in files:
                output += "Results for {}:\n".format(f)
                output += check_pdbqt_order(os.path.join(pdbid, "ledock", f))
                output += "*"*80 + "\n"
    with open("checking_intra.log", "w") as outstream:
        outstream.write(output)

def check_order_inter(database):
    output = ""
    all_substrates = set()
    for uniprot in database["UNIPROT"].values():
        for pdbid in uniprot["pdbs"]:
            for substrate in uniprot["substrates"]:
                all_substrates.add(substrate)
    conformers = _count_conformers("../01_fpocket/{}".format(pdbid))
    for substrate in all_substrates:
        pdbs = []
        for uniprot in database["UNIPROT"].values():
            for pdbid in uniprot["pdbs"]:
                if substrate in uniprot["substrates"]:
                    pdbs.append(pdbid)
        output += "Checking for substrate: {}\n".format(substrate)
        if len(pdbs) > 1:
            sname = "{}_c0_@_{}_LD.pdbqt".format(pdbs[0], substrate)
            order_pdbqt = get_pdbqt_models(os.path.join(pdbs[0], "ledock", sname))[0]
            for pdb in pdbs:
                for i in range(conformers):
                    _name = "{}_c{}_@_{}_LD.pdbqt".format(pdb, i, substrate)
                    to_test = get_pdbqt_models(os.path.join(pdb, "ledock", _name))[0]
                    msg = order_pdbqt.check_atom_order(to_test)
                    if msg is not None:
                        output += msg
                    to_test = None
        output += "*"*80 + "\n"
    with open("checking_inter.log", "w") as outstream:
        outstream.write(output)


if __name__ == "__main__":
    # python3 06_check_conversions.py json_database.json
    with open(sys.argv[1]) as in_stream:
        database = json.load(in_stream)
    check_order_intra(database)
    check_order_inter(database)
