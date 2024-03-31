#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import time
import collections
import json

from libs import paths
from libs import utils
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        substrates_filter_distances()
        logger.info("Execution time: %d s", time.time() - start_time)


def substrates_filter_distances():
    utils.transform_files(
        paths.substrates_clusters(), substrates_cluster_filter,
        TASK["substrates"]["cluster_filter"]["input"],
        TASK["substrates"]["cluster_filter"]["backup"])


def substrates_cluster_filter(input_files, output_files):
    input_file = input_files.pop()
    output_file = output_files.pop()
    logger.info("Filtering clusters: %s -> %s", input_file, output_file)

    fragments = utils.load_json_lines(input_file)

    groups_by_name = collections.defaultdict(list)
    for fragment in fragments:
        groups_by_name[fragment["name"]].append(fragment)

    # Add best molecule for each group using Coulomb value.
    # There are some issues with this approach, since the lower Coulomb value is not
    # always the fragment with shortest distances between atoms.
    filtered_fragments = []
    for group in groups_by_name.values():
        value = min([item["coulomb_score"] for item in group])
        pass_from_group = [item for item in group
                           if item["coulomb_score"] == value]
        filtered_fragments.extend(pass_from_group)
    
    with open(output_file, "w") as out_stream:
        for fragment in filtered_fragments:
            json.dump([fragment], out_stream)
            out_stream.write("\n")


if __name__ == "__main__":
    main()
