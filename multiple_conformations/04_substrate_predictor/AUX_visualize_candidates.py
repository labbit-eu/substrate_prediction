#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  care_visualize_candidates.py
#
#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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
#


import os
import sys
import json
import gzip
import logging
import collections


from libs import utils
from libs import paths
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)


def make_pymol_script(proteins, substrates):
    script = "from pymol import cmd\n"
    for key, value in substrates.items():
        script += "cmd.load('../../candidates/docking/subdir0001/{}')\n".format(key)
        script += "hide sticks\n"
        script += "show lines\n"
        prefix = "S{}_".format(key[:-6])
        script += "split_states {}, prefix={}\n".format(key[:-6], prefix)
        for i in range(1, 601):
            if (i-1) not in value:
                script += "delete {}{:0>4}\n".format(prefix, i)
            else:
                pass
        script += "group {}, {}*\n".format(prefix[:-1], prefix)
        script += "delete {}\n".format(key[:-6])
    for prot in proteins:
        prot_name = prot.split('-')
        prot_name = "{}-{}_{}-{}_{}".format(prot_name[0], prot_name[1],
                                            prot_name[2], prot_name[3],
                                            prot_name[4])
        script += "cmd.load('../../proteins/{}')\n".format(prot_name)
        script += "color green, name Cl\n"
    script += "hide spheres\n"
    script += "hide cartoon\n"
    script += "show ribbon\n"
    return script


def main():
    utils.init_logging()
    args = initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    selected_filters = [int(f) for f in args.filters.split(",")]
    for i in selected_filters:
        cand_files = [f for f in os.listdir(paths.candidates_filtered())
                      if f[:3] == str(i).zfill(3)]
        candidates = []
        for cand in cand_files:
            inpath = os.path.join(paths.candidates_filtered(), cand)
            with gzip.open(inpath, 'rt') as instream:
                for line in instream:
                    candidates.extend(json.loads(line))
        proteins = set()
        mappings = set()
        substrates = collections.defaultdict(list)
        for item in candidates:
            prot = item["ref"]["mapping"][9:] + ".pdbqt"
            proteins.add(prot)
            mappings.add(item["ref"]["mapping"])
            substrates[item["name"]+".pdbqt"].append(item["model_index"])
        out_path = os.path.join(paths.task_root(), "graphs", "candidates")
        os.makedirs(out_path, exist_ok=True)
        out_path = os.path.join(out_path, cand_files[0][:cand_files[0].find(".")]+".pml")
        pymscript = make_pymol_script(proteins, substrates)
        for mapping in mappings:
            map_path = os.path.join(paths.knowledge_base_path(), "mapping",
                                    mapping+".json")
            with open(map_path, 'r') as mapstream:
                jmap = json.load(mapstream)
            sel_name = jmap["ref_id"][9:].split('-')
            sel_name = "{}-{}_{}-{}_{}".format(sel_name[0], sel_name[1],
                                               sel_name[2], sel_name[3],
                                               sel_name[4])
            for res in jmap["properties"]["fasta_to_mcsa"][0]["value"].keys():
                pymscript += "show sticks, {} and resi {}\n".format(sel_name, res)
            pymscript += "color tv_red, name *O*\n"
        pymscript += "delete sele\n"
        with open(out_path, 'w') as outstream:
            outstream.write(pymscript)


if __name__ == "__main__":
    main()
