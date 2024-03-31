#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  care_graph_step_13.py
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
import json
import logging
import collections


from libs import utils
from libs import paths
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)


def make_pymol_script2(proteins, substrates):
    script = "from pymol import cmd\n\n"
    subs = list(substrates.keys())
    subs.sort()
    for key in subs:
        prefix = key[:-6]
        prefix = prefix.replace('(', 'p').replace(')', 'q').replace(',', '_')
        tag = key[:-6].replace('(', 'p').replace(')', 'q').replace(',', '_')
        script += "cmd.load('../../complexes/{}', '{}')\n".format(key, tag)
        script += "cmd.color('gray50', '{} and elem c')\n".format(tag)
        script += "cmd.hide('sticks', '{}')\n".format(tag)
        script += "cmd.show('lines', '{}')\n".format(tag)
        for i in substrates[key]:
            script += "cmd.create('{}_{}', '{}', {}, 1)\n".format(tag, i+1, tag, i+1)
        script += "cmd.group('{}', '{}*')\n".format(prefix.replace('-A', ''), prefix)
        script += "cmd.delete('{}')\n".format(tag)
    proteins = list(proteins)
    proteins.sort()
    for prot in proteins:
        script += "cmd.load('../../proteins/{}')\n".format(prot)
    script += "cmd.hide('spheres', '{}')\n".format(prot)
    script += "cmd.hide('cartoon', '{}')\n".format(prot)
    script += "cmd.show('ribbon', '{}')\n".format(prot)
    return script

def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    clu_files = [f for f in os.listdir(paths.substrates_clusters()) if f[-6:] == ".jsonl"]
    for clu in clu_files:
        inpath = os.path.join(paths.substrates_clusters(), clu)
        cluster = []
        with open(inpath, 'r') as instream:
            for line in instream:
                cluster.extend(json.loads(line))
        proteins = set()
        mappings = set()
        substrates = collections.defaultdict(set)
        for item in cluster:
            substrate = item["ref"]["substrate"][11:]
            prot = item["name"].replace(substrate, "")[:-1] + ".pdbqt"
            proteins.add(prot)
            mappings.add(item["ref"]["mapping"])
            substrates[item["name"]+".pdbqt"].add(item["model_index"])
        out_path = os.path.join(paths.task_root(), "graphs", "clusters")
        os.makedirs(out_path, exist_ok=True)
        out_path = os.path.join(out_path, clu[:-6]+".pml")
        pymscript = make_pymol_script2(proteins, substrates)
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
                pymscript += "select {} and resi {}\n".format(sel_name, res)
                pymscript += "show sticks, sele\n"
            pymscript += "select name *O*\n"
            pymscript += "color tv_red, sele\n"
        pymscript += "delete sele\n"
        with open(out_path, 'w') as outstream:
            outstream.write(pymscript)


if __name__ == "__main__":
    main()
