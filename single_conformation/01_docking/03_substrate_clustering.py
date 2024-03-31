#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  03_substrate_clustering.py
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
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering
from rdkit.Chem.Fingerprints import FingerprintMols
from scipy.cluster.hierarchy import dendrogram


def get_substrate_list(database, subspath):
    substrates = set()
    for uniprot in database["UNIPROT"].values():
        for subs in uniprot["substrates"]:
            temp_sub = "{}/{}.sdf".format(subspath, subs)
            substrates.add(temp_sub)
    return list(substrates)

def plot_dendrogram(model, outfile):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    # Plot the corresponding dendrogram
    plt.figure()
    dendrogram(linkage_matrix, orientation="left", color_threshold=1.32)
    plt.savefig(outfile, dpi=300, format="png")

def make_pymol_script(labels, files):
    subs = [os.path.basename(f) for f in files]
    subs = [s.replace('(', 'p').replace(')', 'q').replace(',', '_')[:-4] for s in subs]
    script = "from pymol import cmd\n"
    for l, f, s in zip(labels, files, subs):
        script += "cmd.load('{}', 'G{}_{}')\n".format(f, l, s)
        script += "cmd.color('gray50', 'G{}_{} and elem c')\n".format(l, s)
    for group in sorted(list(set(labels))):
        script += "cmd.group('cluster_{}', 'G{}_*')\n".format(group, group)
    with open("substrates_clusters.pym", "w") as fout:
        fout.write(script)


def main(datafile):
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    sdf_files = get_substrate_list(database, "substrates")
    sdf_files.sort()
    molecules = [None] * len(sdf_files)
    for i in range(len(sdf_files)):
        molecules[i] = Chem.SDMolSupplier(sdf_files[i])
    fingerprints = [FingerprintMols.FingerprintMol(x[0]) for x in molecules]
    model = AgglomerativeClustering(n_clusters=5, affinity="euclidean", compute_distances=True,
                                   linkage="ward", distance_threshold=None)
    fingerprint_matrix = np.zeros((len(fingerprints), len(fingerprints)))
    for i in range(len(sdf_files)):
        for j in range(len(sdf_files)):
            fingerprint_matrix[i][j] = DataStructs.FingerprintSimilarity(fingerprints[i], fingerprints[j])
    model.fit(fingerprint_matrix)
    np.savetxt("distance_matrix.txt", fingerprint_matrix)
    plot_dendrogram(model, "dendrogram.png")
    groups = [[] for i in range(len(set(model.labels_)))]
    for i in range(len(molecules)):
        groups[model.labels_[i]].append(molecules[i][0])
    tojson = defaultdict(list)
    legend = "ID\tCluster\tFile\n"
    for i, s in enumerate(sdf_files):
        legend += "{}\t{}\t{}\n".format(i, model.labels_[i], os.path.basename(s))
        print(i, s, model.labels_[i])
        tojson[int(model.labels_[i])].append(os.path.basename(sdf_files[i])[:-4])
    with open("substrate_ids.txt", "w") as f:
        f.write(legend)
    make_pymol_script(model.labels_, sdf_files)
    clusters = {}
    for i, group in enumerate(groups):
        common_subs = rdFMCS.FindMCS(group, matchValences=True)
        clusters[i] = (common_subs, group)
        print("Cluster: ", i, common_subs.smartsString)
        smrt_mol = Chem.MolFromSmarts(common_subs.smartsString)
        with open("cluster_{}.sdf".format(i), "w") as outstream:
            outstream.write(Chem.MolToMolBlock(smrt_mol))
            outstream.write("\n> <SMILES>\n")
            outstream.write(Chem.MolToSmiles(smrt_mol))
            outstream.write("\n\n$$$$\n")
    commall = rdFMCS.FindMCS([m[0] for m in molecules], matchValences=True)
    print("Common to all: ", commall.smartsString)
    with open("cluster_all.sdf", "w") as outstream:
        smrt_mol = Chem.MolFromSmarts(commall.smartsString)
        outstream.write(Chem.MolToMolBlock(smrt_mol))
        outstream.write("\n> <SMILES>\n")
        outstream.write(Chem.MolToSmiles(smrt_mol))
        outstream.write("\n\n$$$$\n")
    with open("clusters.json", "w") as outstream:
        json.dump(tojson, outstream, indent=2, sort_keys=True)

# python3 03_substrate_clustering.py json_database.json
main(sys.argv[1])
