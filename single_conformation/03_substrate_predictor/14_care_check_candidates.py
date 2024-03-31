#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  care_check_candidates.py
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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
import time
import json
import logging
import re

from libs import utils
from libs import paths
from libs.configuration import initialize_configuration, TASK

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        check_candidates()
        logger.info("Execution time: %d s", time.time() - start_time)


def check_candidates():
    filter_id = None
    filter_data = []
    visited = set()
    for file_name in sorted(os.listdir(paths.candidates_filtered())):
        name = file_name[:file_name.find("json") - 1]
        if filter_id != name:
            if not _is_filter_empty(filter_id):
                print_details(filter_data, filter_id)
            filter_data = utils.load_json_lines(os.path.join(paths.candidates_filtered(), file_name), gzipped=True)
            filter_id = name
        else:
            filter_data.extend(utils.load_json_lines(os.path.join(paths.candidates_filtered(), file_name), gzipped=True))


def _is_filter_empty(name):
    if name is None:
        return
    with open(os.path.join(paths.filters(), name+".json"), "r") as f:
        in_filter = json.load(f)
    return len(in_filter["filters"]) == 0


def print_details(data, name):
    if len(data) == 0:
        return
    print("-" * 80)
    print("Results for filter id:", name)
    print()
    print("\tNumber of docked conformations:", len(data))
    identified_candidates = set()
    for fragment in data:
        identified_candidates.add(fragment["name"])
    candidates = []
    unidentified_candidates = []
    with open(TASK["candidates"]["prepare"]["manual"], "r") as f:
        for line in f:
            if len(line) > 0:
                candidates.append(line.strip())
    identified = 0
    for c in candidates:
        if c in identified_candidates:
            identified += 1
        else:
            unidentified_candidates.append(c)
    print("\tIdentified candidates:", identified, "/", len(candidates))
    if len(unidentified_candidates) > 0:
        print("\tUndentified candidates:", ", ".join(unidentified_candidates))


if __name__ == "__main__":
    main()
