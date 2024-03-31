#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import json
import shutil


CONF_PATTERN = """[task]
name = {}
knowledge_base.pdb = {}
docking.box = {:4.2f}, {:4.2f}, {:4.2f}
docking.cpu = 10
substrates.cluster_filter.input = ^[0-9]+$
substrates.accessibility.input = ^[0-9]+$
substrates.filter_create.input = ^[0-9]+$
candidates.prepare.manual = {}
candidates.prepare.logP = -10, 8
candidates.prepare.MW = 0, 500
candidates.prepare.TORSDOF = 0, 15
candidates.prepare.max_molecules = 500
candidates.filter_distances.input = ^[0-9]+$
candidates.accessibility.input = ^[0-9]+$
candidates.filter_accessibility.input = ^[0-9]+$
candidates.letters = _A, _B, _C
candidates.docking_mode = manual
substrates.accessibility.save_visualisation = false
substrates.cluster.merge_threshold = 1.0
distances.selected_residues = LIST OF ACTIVE SITE RESIDUES

[env]
vina = PATH_TO_AUTODOCK_VINA
python_mgl = PATH_TO_MGLTOOLS/bin/pythonsh
mgl = PATH_TO_MGLTOOLS/MGLToolsPckgs/AutoDockTools/
babel = PATH_TO_OPENBABEL
python_pymol = PATH_TO_PYMOL/bin/python
clustal_omega = PATH_TO_CLUSTAL_OMEGA
num_parallel_cpu = 1
console_log_level = information
"""

def get_candidates_cids(database):
    allcids = {}
    for p in database["UNIPROT"].values():
        for s,c in p["substrates_cids"].items():
            if s in allcids:
                print("already here", s, c, allcids[s])
                assert c == allcids[s]
            else:
                allcids[s] = c
    return allcids

def make_candidates_list(base_candidates, allcids, data_clusters, cluster, outpath):
    with open(outpath, "w") as outstream:
        for bc in base_candidates:
            outstream.write(str(bc))
            outstream.write("\n")
        to_add = data_clusters[str(cluster)]
        for bc in to_add:
            if bc.endswith("_LD"):
                continue
            outstream.write(str(allcids[bc]))
            outstream.write("\n")

def main():
    candidates_cids = LIST_OF_SUBSTRATE_CIDS_OF_CANDIDATES_TO_TEST
    clusters = 5
    query_protein = "PROTEIN PDB TO QUERY"
    datafile = os.path.join("..", "02_build_dbs", "dehalogenase.json")
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    data_clu_path = os.path.join("..", "02_build_dbs", "clusters.json")
    with open(data_clu_path) as in_stream:
        data_clusters = json.load(in_stream)
    allcids = get_candidates_cids(database)
    box_center = database["DOCKING_BOX"]["size"]
    os.makedirs("config", exist_ok=True)
    
    for clu in range(clusters):
        candidates_path = os.path.join("config", "{}_clu{}_candidates.txt".format(query_protein, clu))
        name = "{}_clu{}".format(query_protein, clu)
        data_from = os.path.join("..", "02_build_dbs", "cluster_{}".format(clu))
        data_to = os.path.join("data", name)
        shutil.copytree(data_from, data_to)
        make_candidates_list(candidates_cids, allcids, data_clusters, clu, candidates_path)
        config_path = os.path.join("config", "{}.ini".format(name))
        config = CONF_PATTERN.format(name, query_protein, box_center[0], box_center[1], box_center[2],
                                     candidates_path)
        with open(config_path, "w") as f:
            f.write(config)


if __name__ == "__main__":
    main()

