#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import time
import os
import shutil

from libs import paths
from libs import filter
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
        candidates_filter_distances()
        logger.info("Execution time: %d s", time.time() - start_time)


def candidates_filter_distances():
    utils.transform_files(
        paths.candidates_filtered(),
        _candidates_filter_distances_file,
        TASK["candidates"]["filter_accessibility"]["input"],
        TASK["candidates"]["filter_accessibility"]["backup"])


def _candidates_filter_distances_file(input_paths, output_paths):
    options = {
        "debug": False
    }

    temp_paths = []
    filter_paths = set()
    for input_path in input_paths:
        name = os.path.basename(input_path)
        temp_path = os.path.join(paths.temp(), name)
        temp_paths.append(temp_path)
        name = name[:name.find(".json")]
        filter_paths.add(os.path.join(paths.filters(), name + ".json"))

    for filter_path in filter_paths:
        filter.filter_file(filter_path, input_paths, temp_paths, options)

    for output_path, temp_path in zip(output_paths, temp_paths):
        shutil.move(temp_path, output_path)


if __name__ == "__main__":
    main()
