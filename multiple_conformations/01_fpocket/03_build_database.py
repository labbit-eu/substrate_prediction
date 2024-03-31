#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import sys
import pickle
import argparse
import numpy as np
from collections import defaultdict
from dataclasses import dataclass, field


@dataclass(order=True, repr=True)
class Pocket:
    pdbfile: str = field(compare=True, hash=True)
    pid: int = field(compare=True, hash=True)
    score: float = field(compare=True, hash=True)
    drug: float = field(compare=True, hash=True)
    residues: set[(str,int)] = field(default_factory=set)
    def __repr__(self):
        return "PDBfile: {}\tID: {}\tScore: {:0>5.4f}\tDrugscore: {:0>5.4f}".format(self.pdbfile, self.pid, self.score, self.drug)

def parse_pockets(infolder):
    entries = [os.path.join(infolder, "pockets", f) for f in os.listdir(os.path.join(infolder, "pockets")) if f.endswith("pdb")]
    entries.sort()
    pockets = []
    _pdbfile = os.path.join(os.path.basename(infolder) + ".pdb")
    for e in entries:
        _residues = set()
        chunks = os.path.basename(e).split("_")
        _pid = int(chunks[0][6:])
        with open(e, "r") as fin:
            for line in fin:
                chunks = line.split()
                if line.startswith("HEADER 0"):
                    _score = float(chunks[-1])
                if line.startswith("HEADER 1"):
                    _drug = float(chunks[-1])
                if line.startswith("ATOM"):
                    _residues.add((chunks[3], int(chunks[4])))
        pockets.append(Pocket(_pdbfile, _pid, _score, _drug, _residues))
    return pockets

def make_database(pockets_path, output_db):
    folders = [f for f in os.listdir(pockets_path)]
    folders.sort()
    pockets_byframe = defaultdict(list)
    for f in folders:
        pockets_byframe[f] = parse_pockets(os.path.join(pockets_path, f))
    with open(output_db, "wb") as fout:
        pickle.dump(pockets_byframe, fout)

def get_best_pockets(pockets):
    best_pockets = []
    used_resids = set()
    for p in pockets:
        res = p.pdbfile.split("_")[1].split(".")[0]
        if res in used_resids:
            continue
        else:
            used_resids.add(res)
            best_pockets.append(p)
    return best_pockets

def print_best_pockets(datafile):
    with open(datafile, "rb") as fin:
        data = pickle.load(fin)
    pockets = []
    for key in data.keys():
        for pkt in data[key]:
            if (("GLU",45) in pkt.residues) and (("LYS",132) in pkt.residues) and \
               (("CYS",166) in pkt.residues):
                pockets.append(pkt)
    pockets.sort(key=lambda p: p.drug*p.score, reverse=True)
    pockets = get_best_pockets(pockets)
    script = "from pymol import cmd\n"
    pdbid = pockets[0].pdbfile[:6]
    script += "cmd.load('../../02_lrip_sampling/{}/{}.tleap.pdb', 'candidate_0')\n".format(pdbid, pdbid)
    for i, p in enumerate(pockets):
        print(p)
        script += "cmd.load('pdbs/{}', 'candidate_{}')\n".format(p.pdbfile.replace("_out", ""), i+1)
    script += "cmd.hide('cartoon', 'all')\n"
    script += "cmd.show('ribbon', 'all')\n"
    script += "cmd.show('sticks', 'resi 45+132+166')\n"
    script += "cmd.zoom('resi 45+132+166')\n"
    with open("proteins.pym", "w") as fout:
        fout.write(script)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        # python3 03_build_database.py fpockets pockets3.1.dat
        make_database(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        # python3 03_build_database.py pockets1.dat
        print_best_pockets(sys.argv[1])
