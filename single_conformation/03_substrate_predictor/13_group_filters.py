#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import time
import json
import logging

from libs import paths
from libs import utils
from collections import defaultdict
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        groups = group_filters()
        filter_group_by_residues(groups)
        logger.info("Execution time: %d s", time.time() - start_time)


def group_filters():
    _filters = os.listdir(paths.filters())
    filters = []
    for f in _filters:
        name = f[:f.find(".json")]
        if re.match(TASK["substrates"]["group_filters"]["input"], name):
            filters.append(f)
    filters.sort(reverse=True)
    f_groups = [[filters.pop()]]
    while filters:
        f = filters.pop()
        if not _similar_to_any(f, f_groups):
            f_groups.append([f])
    return f_groups


def _similar_to_any(input_filter, groups):
    for group in groups:
        if _has_same_parameters(input_filter, group[0]):
            group.append(input_filter)
            return True
    return False


def _has_same_parameters(input_filter, group):
    with open(os.path.join(paths.filters(), input_filter), "r") as f:
        in_filter = json.load(f)
    with open(os.path.join(paths.filters(), group), "r") as f:
        to_compare = json.load(f)
    if len(in_filter["filters"]) != len(to_compare["filters"]):
        return False
    for in_fil in in_filter["filters"]:
        flag = False
        for out_fil in to_compare["filters"]:
            if in_fil["keys"][0] == out_fil["keys"][0]:
                flag = True
        if not flag:
            return False
    return True


def filter_group_filters(groups):
    num_keys = len(TASK["distances"]["selected"]) * 2
    new_groups = defaultdict(list)
    for i, group in enumerate(groups):
        for element in group:
            with open(os.path.join(paths.filters(), element), "r") as f:
                in_filter = json.load(f)
            if len(in_filter["filters"]) != num_keys:
                continue
            for in_fil in in_filter["filters"]:
                if in_fil["prop"] == "accessibility":
                    if in_fil["intervals"][0][0] < 0.5:
                        break
                elif in_fil["prop"] == "distances":
                    if in_fil["intervals"][0][1] > 5.0:
                        break
            else:
                new_groups[i].append(element)
    for gid, group in new_groups.items():
        print("The most likely filter to be correct is defined by:", group)


def _all_residues_present(input_filter):
    with open(os.path.join(paths.filters(), input_filter), "r") as f:
        in_filter = json.load(f)
    flag = True
    for residue in TASK["distances"]["selected_residues"]:
        for in_fil in in_filter["filters"]:
            if residue in in_fil["keys"][0]:
                break
        else:
            flag = False
    return flag


def filter_group_by_residues(groups, unknown_roles=False):
    num_keys = TASK["distances"]["selected_residues"]
    new_groups = defaultdict(list)
    if unknown_roles:
        thresholds = (0.3, 6.0)
    else:
        thresholds = (0.5, 5.0)
    for i, group in enumerate(groups):
        if not _all_residues_present(group[0]):
            continue
        for element in group:
            with open(os.path.join(paths.filters(), element), "r") as f:
                in_filter = json.load(f)
            for in_fil in in_filter["filters"]:
                if in_fil["prop"] == "accessibility":
                    if in_fil["intervals"][0][0] < thresholds[0]:
                        break
                elif in_fil["prop"] == "distances":
                    if in_fil["intervals"][0][1] > thresholds[1]:
                        break
            else:
                new_groups[i].append(element)
    for gid, group in new_groups.items():
        print("The most likely filter to be correct is defined by:", group)


if __name__ == "__main__":
    main()
