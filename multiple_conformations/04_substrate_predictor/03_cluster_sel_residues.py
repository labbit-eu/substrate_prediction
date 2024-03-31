#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import json
import time
import hdbscan
import logging
import numpy as np
import sklearn.cluster as skcluster

from libs import utils
from libs import paths
from collections import defaultdict
import matplotlib.pyplot as plt
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)
logging.getLogger("matplotlib.font_manager").disabled = True


def main():
    utils.init_logging()
    initialize_configuration()
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    distances = TASK["distances"]["selected"]
    with utils.FileLogging(paths.logs(), __file__):
        substrates_cluster()
        logger.info("Execution time: %d s", time.time() - start_time)


def substrates_cluster():
    in_path = paths.substrates_eval()
    out_dir = paths.substrates_clusters()
    logger.info("Creating files in `%s` from file `%s`", out_dir, in_path)
    create_clusters_for_file(in_path, out_dir)


def create_clusters_for_file(input_file, output_dir):
    fragments = utils.load_json_lines(input_file)
    clusters = _create_clusters(fragments)
    clusters.sort(key=lambda cluster: len(cluster), reverse=True)
    _save_clusters(clusters, output_dir)


def _create_clusters(fragments):
    options = TASK["substrates"]["cluster"]
    clusters = _cluster_by_distances_rmsd(fragments, options["merge_threshold"], "kmeans")
    return clusters


# region clustering

def _cluster_by_distances_rmsd(fragments, merge_threshold, cluster_method="hdbscan"):
    logger.info("Feeding HDBSCAN with computed distances (%d fragments)", len(fragments))
    _remove_unused_residues(fragments)
    distances = _create_distance_descriptors(fragments)
    hbdclusters = hdbscan.HDBSCAN(min_cluster_size=10, cluster_selection_epsilon=merge_threshold).fit(distances)
    
    if cluster_method == "hdbscan":
        logger.info("HDBSCAN clustering ...")
        merge_records = hbdclusters.condensed_tree_.to_numpy()
        clusters = _calculate_fragment_clusters_hdbscan(merge_records, fragments, merge_threshold)
    else:
        logger.info("KMeans clustering ...")
        hbdclusters = len(set(hbdclusters.labels_))
        clusterer = skcluster.KMeans(n_clusters=hbdclusters).fit(distances)
        clusters = _no_join_clusters(clusterer.labels_, fragments)
    
    logger.info("Setting clusters")
    return clusters


def _no_join_clusters(labels, fragments):
    hdb_clusters = [[] for l in set(labels)]
    for l, fragment in zip(labels, fragments):
        hdb_clusters[l].append(fragment)
    return hdb_clusters


def _remove_unused_distances(fragments):
    for fragment in fragments:
        new_distances = []
        for distance in fragment["distances"]:
            if _is_selected(distance["name"]):
                new_distances.append(distance)
        fragment["distances"] = new_distances


def _remove_unused_residues(fragments):
    for fragment in fragments:
        new_distances = []
        for distance in fragment["distances"]:
            if _is_residue_selected(distance["name"]):
                new_distances.append(distance)
        fragment["distances"] = new_distances


def _create_distance_descriptors(fragments):
    """
    Use first fragment to define order of distances, then order distances
    for all fragment according to the ordering.

    The output are new ordered distances for each fragment.
    frag_1 = [1.8, 0.1, 3.1, 4.4]
    dists = [frag_1, frag_2, frag_3, frag_4]
    """
    output = []
    keys = [item["name"] for item in fragments[0]["distances"]]
    for fragment in fragments:
        distances = [0] * len(keys)
        for item in fragment["distances"]:
            index = keys.index(item["name"])
            if index == -1:
                raise RuntimeError("Unexpected distance keys in: {}".format(fragment["name"]))
            distances[index] = item["val"]
        output.append(distances)
    return output


def _is_selected(name):
    for sel_dist in TASK["distances"]["selected"]:
        if sel_dist[0] == name[1] and sel_dist[1] == name[2] and \
           sel_dist[2] == name[3] and sel_dist[3] == name[4]:
            return True
    return False


def _is_residue_selected(name):
    if name[1] in TASK["distances"]["selected_residues"]:
        return True
    return False


def _calculate_fragment_clusters_hdbscan(hdb_results, fragments, merge_threshold):
    hdb_clusters = [cluster for cluster in hdb_results]
    hdb_clusters.sort(key=lambda x: x[2], reverse=True)
    clusters = [[fragment] for fragment in fragments]
    pos = 0
    for cluster in hdb_clusters:
        if cluster[2] < merge_threshold:
            break
        pos += 1
    for i in range(len(clusters), len(hdb_clusters)+1):
        clusters.append([])
    clustering = hdb_clusters[:pos]
    clustering.sort(key=lambda x: x[0], reverse=True)
    for record in clustering:
        left_index = int(record[0])
        right_index = int(record[1])
        clusters[left_index] += clusters[right_index]
        clusters[right_index] = None
    # Filter out None ie merged clusters.
    clusters = [cluster for cluster in clusters if cluster is not None]
    clusters = [cluster for cluster in clusters if len(cluster)>0]
    return clusters


# endregion clustering


def _save_clusters(clusters, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    counter = 0
    for fragments in clusters:
        if len(fragments) < TASK["substrates"]["cluster"]["min_cluster_size"]:
            continue
        out_path = os.path.join(output_dir, str(counter).zfill(3) + ".jsonl")
        with open(out_path, "w") as out_stream:
            for fragment in fragments:
                json.dump([fragment], out_stream)
                out_stream.write("\n")
        counter += 1


if __name__ == "__main__":
    main()

