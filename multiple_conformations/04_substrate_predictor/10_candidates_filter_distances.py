#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import time
import os

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
    logger.info("Filtering candidates: %s -> %s",
                str(paths.candidates_eval()),
                str(paths.candidates_filtered()))

    utils.transform_files(
        paths.filters(),
        _candidates_filter_distances_for_filter,
        TASK["candidates"]["filter_distances"]["input"],
        None)


def _candidates_filter_distances_for_filter(filter_paths, _):
    options = {
        "ignore-props": ["accessibility"],
        "debug": False
    }
    in_paths = paths.candidates_eval()

    for filter_path in filter_paths:
        out_paths = []
        for in_path in in_paths:
            out_name = os.path.basename(filter_path) + "_{}.gz".format(paths.get_candidates_eval_id(in_path))
            out_paths.append(os.path.join(paths.candidates_filtered(), out_name))

        filter.filter_file(filter_path, in_paths, out_paths, options)


if __name__ == "__main__":
    main()
