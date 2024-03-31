#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os
import gzip
import multiprocessing as mp

# from candidates_expected import VERIFIED, EXPECTED
from libs.configuration import ENV, _SPECIAL_RESIDUES

logger = logging.getLogger(__name__)


def filter_file(filter_path, input_paths, output_paths, options):
    filters, used_props = _load_filters(filter_path, options)
    input_size = 0
    output_size = 0

    pool = mp.Pool(processes=ENV["num_parallel_cpu"])
    processing = []

    for input_path, output_path in zip(input_paths, output_paths):
        processing.append(pool.apply_async(_filter_file, args=(input_path, output_path, options, filters, used_props)))

    for p in processing:
        tmp = p.get()
        input_size += tmp[0]
        output_size += tmp[1]

    if options.get("debug", True):
        _print_filters_definition(filters)

    logger.info("Filtering for %s done %d -> %d", os.path.basename(filter_path), input_size, output_size)


def _filter_file(input_path, output_path, options, filters, used_props):
    input_size = 0
    output_size = 0

    with gzip.open(output_path, "wt") as out_stream:
        for batch in _iter_fragment_batches(input_path):
            filtered_fragments = [
                fragment for fragment in batch
                if _do_pass_filters(filters, used_props, fragment, options)]
            input_size += len(batch)
            output_size += len(filtered_fragments)
            if len(filtered_fragments) == 0:
                continue
            logger.debug("Filtering batch %s: %d -> %d", batch[0]["name"],
                         len(batch), len(filtered_fragments))
            json.dump(filtered_fragments, out_stream)
            out_stream.write("\n")

    return [input_size, output_size]


def _load_filters(filter_path, options):
    with open(filter_path) as in_stream:
        filters = json.load(in_stream)["filters"]

    # Ignore filters.
    if "ignore-props" in options:
        filters = [f for f in filters
                   if f["prop"] not in options["ignore-props"]]

    # Apply tolerance to intervals.
    for filter_ in filters:
        if "tolerance" not in filter_:
            continue
        tolerance = filter_["tolerance"]
        for index, interval in enumerate(filter_["intervals"]):
            if filter_["prop"] == "distances":
                new_interval = [interval[0], interval[1] + tolerance]
            elif filter_["prop"] == "accessibility":
                new_interval = [max(interval[0] - tolerance, 0), interval]
            else:
                raise RuntimeError("Unknown property" + filter_["prop"])
            filter_["intervals"][index] = new_interval

    # JSON can not store tuples, only arrays. So we need to convert this back.
    for filter_ in filters:
        filter_["keys"] = list(map(tuple, filter_["keys"]))

    used_props = set([f["prop"] for f in filters])
    return filters, used_props


def _iter_fragment_batches(input_path):
    with gzip.open(input_path) as in_stream:
        for line in in_stream:
            yield json.loads(line.decode('utf-8'))


def _do_pass_filters(filters, used_props, fragment, options):
    try:
        values = {
            prop: _array_to_directory(fragment, prop)
            for prop in used_props
        }
    except RuntimeError:
        logger.exception("Allowing invalid molecule to pass.")
        return True
    for filter_ in filters:
        try:
            
            outcome, _ = eval_filter_from_definition(
                filter_, values["distances"], values[filter_["prop"]])

            if not outcome and options.get("debug", False):
                _print_debug(fragment, filters, values)

        except KeyError:
            logger.exception("Fragment: " + str(fragment))

            return True
        if not outcome:
            return False
    return True


def _print_debug(fragment, filters, values):
    label = None
    if fragment["name"] in VERIFIED:
        label = "* {}, {}".format(
            fragment["name"], str(fragment["model_index"]))
    elif fragment["name"] in EXPECTED:
        label = "+ {}, {}".format(
            fragment["name"], str(fragment["model_index"]))
    else:
        label = "  {}, {}".format(
            fragment["name"], str(fragment["model_index"]))

    if label is not None:
        eval_debug_filters_from_definitions(
            filters, values, label)


def _array_to_directory(fragment, prop_name):
    if prop_name not in fragment:
        logger.error("Missing property '%s' on : %s", prop_name, str(fragment))
        raise RuntimeError("Can't prepare values for fragment.")
    directory = {}
    for item in fragment[prop_name]:
        directory[tuple(item["name"])] = item["val"]
    return directory


# region filters evaluation

def eval_filter_from_definition(filter_, distances, values):
    if filter_["@type"] == "ordered":
        return _eval_ordered(filter_, values)
    elif filter_["@type"] == "single":
        return _eval_single(filter_, distances, values)
    else:
        message = "Unsupported filter: " + filter_["@type"] + \
                  " : " + filter_["prop"]
        raise RuntimeError(message)


def _eval_ordered(definition, values):
    values = sorted([values[key] for key in definition["keys"]])
    outcomes = []
    errors = []
    for index, interval in enumerate(definition["intervals"]):
        value = values[index]
        outcome, error = _check_value_in_interval(interval, value)
        outcomes.append(outcome)
        errors.append(error)
    return min(outcomes), errors


def _check_value_in_interval(interval, value):
    if interval[0] > value:
        return False, interval[0] - value
    if interval[1] < value:
        return False, value - interval[1]
    return True, 0


def _eval_single(filter_, distances, values):
    if filter_["prop"] == "distances":
        return _eval_single_distances(filter_, values)
    elif filter_["prop"] == "accessibility":
        return _eval_single_accessibility(filter_, distances, values)


def _eval_single_distances(definition, values):
    value = min([values[key] for key in definition["keys"]])
    interval = definition["intervals"][0]
    outcome, error = _check_value_in_interval(interval, value)
    return outcome, [error]


def _eval_single_accessibility(definition, distances, values):
    value = select_accessibility_to_use(definition, distances, values)
    interval = definition["intervals"][0]
    outcome, error = _check_value_in_interval(interval, value)
    return outcome, [error]


def select_accessibility_to_use(definition, distances, values):
    # Select atom with minimum distance.
    min_distance = min([distances[key] for key in definition["keys"]])
    # print("FILTER.PY1:", distances)
    # print("FILTER.PY2:", definition["keys"], min_distance)
    # Select maximum from values that have minimum distance.
    value = max([values[key] for key in definition["keys"]
                 if distances[key] == min_distance])
    return value


# endregion

# region filters evaluation with debug output

def eval_debug_filters_from_definitions(filters, values_dict, label):
    all_outcomes = []
    all_errors = []
    pass_all_filters = True

    for filter_ in filters:
        values = values_dict[filter_["prop"]]
        distances = values_dict["distances"]
        outcome, errors = eval_filter_from_definition(filter_, distances, values)
        all_outcomes.append(outcome)
        all_errors.append(errors)
        pass_all_filters &= outcome

    if pass_all_filters:
        return

    print(" ", label)
    for filter_, outcome, error in zip(filters, all_outcomes, all_errors):
        if outcome:
            continue

        if filter_["@type"] == "single" and filter_["prop"] == "accessibility":
            _print_debug_single_accessibility(filter_, values_dict)
        else:
            _print_debug_default(filter_, values_dict)


def _print_debug_default(filter_, values_dict):
    values = [str(values_dict[filter_["prop"]][key])
              for key in filter_["keys"]]

    print("      actual: {:<30}".format(", ".join(values)),
          " expected: {:<30}".format(", ".join(map(str, filter_["intervals"]))),
          "{:<15}".format(filter_["prop"]),
          ", ".join(["{} {} {} {}".format(k[0], k[1], k[3], k[4])
                     for k in filter_["keys"]]))


def _print_debug_single_accessibility(filter_, values_dict):
    values = [
        "{} ({})".format(values_dict[filter_["prop"]][key],
                         values_dict["distances"][key])
        for key in filter_["keys"]]

    print("      actual: {:<30}".format(", ".join(values)),
          " expected: {:<30}".format(", ".join(map(str, filter_["intervals"]))),
          "{:<15}".format(filter_["prop"]),
          ", ".join(["{} {} {} {}".format(k[0], k[1], k[3], k[4])
                     for k in filter_["keys"]]))


def _print_filters_definition(filters):
    print()
    print("Filters:")
    for filter_ in filters:
        print(" ", filter_["@type"], filter_["prop"])
        print("   keys:     ", ", ".join(["{} {} {} {}".format(
            k[0], k[1], k[3], k[4]) for k in filter_["keys"]]))
        print("   intervals: ", ", ".join(["[{},{}]".format(
            interval[0], interval[1]) for interval in filter_["intervals"]]))
    print()

# endregion
