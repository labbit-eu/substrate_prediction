#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  01_build_dbs.py
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
from collections import defaultdict


def get_substrate_list(database, subspath):
    substrates = set()
    for uniprot in database["UNIPROT"].values():
        for subs in uniprot["substrates"]:
            temp_sub = "{}/{}.sdf".format(subspath, subs)
            substrates.add(temp_sub)
    return list(substrates)

def build_subdatabases(database, clusters, conformers=6):
    for clu, cluster in clusters.items():
        clu_path = "cluster_{}".format(clu)
        conf_name = database["QUERY_PROTEIN"]["PDB"]
        os.makedirs(clu_path, exist_ok=True)
        _query_prot = {"NAME":conf_name, "PDB":conf_name, "substrates":database["QUERY_PROTEIN"]["substrates"],
                       "substrates_cids":database["QUERY_PROTEIN"]["substrates_cids"]}
        _query_prot["active_sites"] = {conf_name:database["QUERY_PROTEIN"]["active_sites"][database["QUERY_PROTEIN"]["PDB"]]}
        _db = {"REF_PDB":database["REF_PDB"], "DOCKING_BOX":database["DOCKING_BOX"],
               "MAPPING_INFO":database["MAPPING_INFO"], "QUERY_PROTEIN":_query_prot}
        _db["UNIPROT"] = defaultdict(dict)
        for u_id, entry in database["UNIPROT"].items():
            for i in range(conformers):
                new_uid = "{}_c{}".format(u_id, i)
                subset = []
                for sub in entry["substrates"]:
                    if sub in cluster:
                        continue
                    subset.append(sub)
                    subset.append(sub+"_LD")
                if len(subset) == 0:
                    continue
                _db["UNIPROT"][new_uid]["substrates"] = subset
                _db["UNIPROT"][new_uid]["pdbs"] = ["{}_c{}".format(entry["pdbs"][0], i)]
                _active_site = {}
                _active_site["{}_c{}".format(entry["pdbs"][0], i)] = entry["active_sites"][entry["pdbs"][0]]
                _db["UNIPROT"][new_uid]["active_sites"] = _active_site
        with open(os.path.join(clu_path, "subdatabase.json"), "w") as dbout:
            json.dump(_db, dbout, indent=2)
        _db.clear()

def create_links(database):
    os.makedirs("docking", exist_ok=True)
    os.makedirs("receptors", exist_ok=True)
    os.makedirs("substrates", exist_ok=True)
    os.makedirs(os.path.join("docking", "logs"), exist_ok=True)
    os.makedirs(os.path.join("docking", "outputs"), exist_ok=True)
    for uniprot in database["UNIPROT"].values():
        # AutoDock vina results
        dock_folder = os.path.join("..", "01_docking", uniprot["pdbs"][0], "outputs")
        dockresults = [f for f in os.listdir(dock_folder)]
        for dr in dockresults:
            os.symlink(os.path.join("..", "..", dock_folder, dr), os.path.join("docking", "outputs", dr))
        log_folder = os.path.join("..", "01_docking", uniprot["pdbs"][0], "logs")
        logresults = [f for f in os.listdir(log_folder)]
        for lr in logresults:
            os.symlink(os.path.join("..", "..", log_folder, lr), os.path.join("docking", "logs", lr))
        # LeDock results
        dock_folder = os.path.join("..", "01_docking", uniprot["pdbs"][0], "ledock")
        dockresults = [f for f in os.listdir(dock_folder) if f.endswith("pdbqt")]
        for dr in dockresults:
            os.symlink(os.path.join("..", "..", dock_folder, dr), os.path.join("docking", "outputs", dr))
    rec_folder = os.path.join("..", "01_docking", "receptors")
    receptors = [f for f in os.listdir(rec_folder) if not f.endswith("_LD.pdb")]
    sub_folder = os.path.join("..", "01_docking", "substrates")
    substrates = [f for f in os.listdir(sub_folder) if f.endswith("sdf") or f.endswith("pdbqt")]
    for r in receptors:
        os.symlink(os.path.join("..", rec_folder, r), os.path.join("receptors", r))
    for s in substrates:
        os.symlink(os.path.join("..", sub_folder, s), os.path.join("substrates", s))
        # LeDock results
        if s.endswith("sdf"):
            os.symlink(os.path.join("..", sub_folder, s), os.path.join("substrates", s.replace(".sdf", "_LD.sdf")))
        else:
            os.symlink(os.path.join("..", sub_folder, s), os.path.join("substrates", s.replace(".pdbqt", "_LD.pdbqt")))

def main(datafile, clustersfile):
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    with open(clustersfile) as in_stream:
        clusters = json.load(in_stream)
    for uniprot in database["UNIPROT"].values():
        for pdbid in uniprot["pdbs"]:
            conformers = 1
    build_subdatabases(database, clusters, conformers)
    create_links(database)

# python3.9 01_build_dbs.py json_database.json clusters.json
main(sys.argv[1], sys.argv[2])
