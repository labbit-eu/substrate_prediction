#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import logging
import time
import os
import collections
import functools
import re
import copy

from libs import file_formats, filter as filter_module
from libs import paths
from libs.data_model import ReadOnlyDatabase, McsaEntry, MappingEntry
from libs import active_site as active_site_lib
from libs import utils
from libs.configuration import initialize_configuration, TASK, ENV
from libs import interactions

logger = logging.getLogger(__name__)


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    start_time = time.time()
    with utils.FileLogging(paths.logs(), __file__):
        substrates_filter_create()
        logger.info("Execution time: %d s", time.time() - start_time)


def substrates_filter_create():
    cluster_dir = paths.substrates_clusters()
    logger.info("Creating filters in '%s' based on data in '%s'",
                paths.filters(), cluster_dir)
    for file_name in sorted(os.listdir(cluster_dir)):
        name = file_name[:file_name.find("jsonl") - 1]
        if not re.match(TASK["substrates"]["filter_create"]["input"], name):
            continue

        in_path = os.path.join(cluster_dir, file_name)
        out_path = os.path.join(paths.filters(), name + ".json")

        create_filters(in_path, out_path, {
            "relax_threshold": TASK["substrates"]["filter_create"]["relax_threshold"],
            "distance_type": TASK["substrates"]["filter_create"]["distance_type"],
            "accessibility_type": TASK["substrates"]["filter_create"]["accessibility_type"],
            "normalization": TASK["substrates"]["filter_create"]["normalization"]
        })


def create_filters(input_path, output_path, options):
    print("-" * 80)
    print("" * 10, input_path, " --> ", output_path)
    logger.info("Creating filters for: %s", os.path.basename(input_path))
    database = ReadOnlyDatabase(paths.database_directory())
    fragments = utils.load_json_lines(input_path)
    filters, metadata = _create_filters(database, fragments, options)
    _save_filters(output_path, {
        "filters": filters,
        "metadata": metadata
    })


def _save_filters(output_path, filters):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out_stream:
        json.dump(filters, out_stream, indent=2)


def _create_filters(database, fragments, options):
    filters = []

    # TODO Reuse output from distances also for accessibility.
    filters.extend(_create_distances_filter_template(
        database, fragments, options))

    _array_to_directory(fragments, "distances")
    _array_to_directory(fragments, "accessibility")
    filters.extend(create_accessibility_filter_template(
        fragments, options))
    filters = _check_filters_dist_access(filters)

    if len(filters) == 0:
        logger.info("No filter templates found!")
        return [], {}

    initialize_filters(filters)
    filters, metadata = relax_filters(fragments, filters, options)

    for filter_ in filters:
        del filter_["$"]

    return filters, {"relaxation": metadata}


def _array_to_directory(fragments, prop_name):
    for fragment in fragments:
        if prop_name not in fragment:
            continue
        directory = {}
        for item in fragment[prop_name]:
            directory[tuple(item["name"])] = item["val"]
        fragment[prop_name] = directory


def _collect_property_values(fragments, prop_name):
    output = collections.defaultdict(list)
    for fragment in fragments:
        for key, value in fragment[prop_name].items():
            output[key].append(value)
    return output


# region create distance filters template


def _create_distances_filter_template(database, fragments, options):
    distances = _select_relevant_distances(database, fragments)
    relevant_distances = interactions.filter_distances_by_type(distances)

    filter_type = options.get("distance_type", "ordered")
    templates = {}
    for item in relevant_distances.keys():
        res_seq, res_name, res_type, res_atom, frag_atom = item
        key = (res_seq, res_name, res_type, res_atom, frag_atom)
        if key not in templates:
            tolerance = active_site_lib.activity_tolerance(res_name, res_type)
            threshold = active_site_lib.activity_threshold(res_name, res_type)
            templates[key] = {
                "@type": filter_type,
                "prop": "distances",
                "keys": [],
                "tolerance": tolerance,
                "$": {
                    "threshold": [0, threshold],
                }
            }
        templates[key]["keys"].append(item)

    return list(templates.values())


def _select_relevant_distances(database, fragments):
    """
    Select distances that should be used in the filter.
    """
    distances = collections.defaultdict(list)
    for fragment in fragments:
        molecule = load_molecule(fragment["name"])[fragment["model_index"]]
        protein = load_protein_atoms(database, fragment["ref"]["mapping"])
        new_distances = interactions.select_reactive_distances(
            fragment, molecule, protein)
        interactions.add_distances(new_distances, distances)
    return distances


# endregion

# region pdbqt loading

@functools.lru_cache(maxsize=16)
def load_protein_atoms(database, mapping_ref):
    """
    :return Atoms that are part of active site and are specified in MCSA record.
    """
    mapping = database.read(MappingEntry, mapping_ref)
    mcsa = database.read(McsaEntry, mapping.mcsa().first().value())

    protein_path = paths.protein_pdbqt_path(mapping, mcsa)
    protein_pdbqt = file_formats.load_pdbqt_file(protein_path)

    active_site = active_site_lib.active_site_from_mcsa(mcsa)
    atoms_by_names = active_site_lib.reacting_atoms(
        protein_pdbqt, mapping, active_site)

    return {
        "reactive_by_name": atoms_by_names
    }


@functools.lru_cache(maxsize=64)
def load_molecule(name):
    path = paths.complexes_pdbqt_path(name)
    return file_formats.load_pdbqt_models_file(path)


# endregion

# region create accessibility filters template

def create_accessibility_filter_template(fragments, options):
    names = set()
    names_to_eval = collections.defaultdict(int)
    for fragment in fragments:
        names.update(fragment["accessibility"].keys())
        for key, value in fragment["accessibility"].items():
            threshold = active_site_lib.activity_threshold(key[1], key[2], "accessibility")
            if value >= threshold:
                names_to_eval[key] += 1
    for key, value in names_to_eval.items():
        if value < len(fragments)*options["relax_threshold"]:
            names.remove(key)
    
    templates = {}
    filter_type = options.get("accessibility_type", "ordered")
    for name in names:
        res_seq, res_name, res_type, res_atom, frag_atom = name
        key = (res_seq, res_name, res_type, res_atom, frag_atom)
        if key not in templates:
            templates[key] = {
                "@type": filter_type,
                "prop": "accessibility",
                "keys": [],
                "$": {
                    "threshold": [0, 1],
                }
            }
        templates[key]["keys"].append(name)

    return list(templates.values())


# endregion

def initialize_filters(filters):
    # TODO We should use min/max values from filters.
    for filter_ in filters:
        if filter_["prop"] == "accessibility":
            start_value = 1
        else:
            start_value = 0
        if filter_["@type"] == "single":
            filter_["intervals"] = [[start_value, start_value]]
            continue
        filter_["intervals"] = [[start_value, start_value]
                               for _ in filter_["keys"]]


# region distance filters relaxation

def relax_filters(fragments, filters, options):
    filters, metadata, to_discard, output = _test_relax_filters(fragments, filters, options)
    while to_discard is not None:
        filters = _remove_filter(filters, to_discard)
        filters, metadata, to_discard, output = _test_relax_filters(fragments, filters, options)
    return filters, metadata


def add_normalize_functions_to_filters(filters, fragments, options):
    values_by_prop = {
        "distances": _collect_property_values(fragments, "distances"),
        "accessibility": _collect_property_values(fragments, "accessibility")
    }

    value_selectors_by_prop = {
        "distances": min,
        "accessibility": max
    }
    for filter_ in filters:
        keys_values = {}
        for key in filter_["keys"]:
            keys_values[key] = values_by_prop[filter_["prop"]][key]
        if filter_["@type"] == "single":
            values = []
            selector = value_selectors_by_prop[filter_["prop"]]
            for index in range(0, len(fragments)):
                value = selector([val[index] for val in keys_values.values()])
                values.append(value)
            filter_["normalize"] = [
                create_normalize_function(values, options)
            ]
        else:
            ordered_values = [[] for _ in filter_["keys"]]
            for index in range(0, len(fragments)):
                values = sorted([val[index] for val in keys_values.values()])
                for index in range(0, len(values)):
                    ordered_values[index].append(values[index])
            filter_["normalize"] = [
                create_normalize_function(values, options)
                for values in ordered_values
            ]


def create_normalize_function(values, options):
    if options["normalization"] == "none":
        return lambda value: value
    elif options["normalization"] == "range":
        rang = max(values) - min(values)
        if rang == 0:
            return lambda value: value
        else:
            return lambda value: value / rang

    message = "Unsupported normalization: " + options["normalization"]
    raise RuntimeError(message)


def relax_filter(filter_, fragment):
    values = fragment[filter_["prop"]]
    if filter_["@type"] == "ordered":
        return relax_ordered(filter_, values)

    if filter_["@type"] == "single":
        if filter_["prop"] == "distances":
            return relax_single_distances(filter_, values)
        elif filter_["prop"] == "accessibility":
            return relax_single_accessibility(
                filter_, fragment["distances"], values)

    message = "Unsupported filter: " + filter_["@type"] + " : " + filter_["prop"]
    raise RuntimeError(message)


def relax_ordered(definition, distances):
    values = sorted([distances[key] for key in definition["keys"]])
    errors = []
    for index, interval in enumerate(definition["intervals"]):
        value = values[index]
        _, error = check_value_in_interval(interval, value)
        errors.append(error)
    non_zero_min = min([e for e in errors if e > 0])
    index_to_relax = errors.index(non_zero_min)

    definition["intervals"][index_to_relax] = relax_interval_max(
        definition["intervals"][index_to_relax], values[index_to_relax])


def relax_interval_max(interval, value):
    return [interval[0], max(interval[1], value)]


def relax_single_distances(definition, values):
    value = min([values[key] for key in definition["keys"]])

    definition["intervals"][0] = relax_interval_max(
        definition["intervals"][0], value)


def relax_single_accessibility(definition, distances, values):
    value = filter_module.select_accessibility_to_use(
        definition, distances, values)

    definition["intervals"][0] = relax_interval_min(
        definition["intervals"][0], value)


def relax_interval_min(interval, value):
    return [min(interval[0], value), interval[1]]


def check_value_in_interval(interval, value):
    if interval[0] > value:
        return False, interval[0] - value
    if interval[1] < value:
        return False, value - interval[1]
    return True, 0


def filter_bad_fragments(fragments, filters):
    """
    We need to get rid of molecules that are over thresholds.
    For this reason we set intervals to thresholds and just evaluate
    the filters.
    """
    threshold_filters = []
    for filter_ in filters:
        threshold_filters.append({
            "@type": filter_["@type"],
            "keys": filter_["keys"],
            "prop": filter_["prop"],
            "intervals": [filter_["$"]["threshold"] for _ in filter_["intervals"]]
        })

    good_fragments = []
    bad_fragments = []
    for fragment in fragments:
        for filter_ in threshold_filters:
            distances = fragment["distances"]
            values = fragment[filter_["prop"]]
            is_ok, _ = filter_module.eval_filter_from_definition(
                filter_, distances, values)
            if not is_ok:
                bad_fragments.append(fragment)
                break
        else:
            good_fragments.append(fragment)

    return good_fragments, bad_fragments


def _test_relax_filters(fragments, filters, options):
    output = ""
    orig_filters = copy.deepcopy(filters)
    add_normalize_functions_to_filters(filters, fragments, options)

    output += "Relaxation on {} input molecules:\n".format(len(fragments))
    pass_threshold = len(fragments) * options["relax_threshold"]
    fragments_to_relax, bad_fragments = filter_bad_fragments(fragments, filters)
    output += "    {} molecules removed from relaxation\n".format(len(bad_fragments))
    output += "    {} molecules left for relaxation\n\n".format(len(fragments_to_relax))

    iteration_counter = 0
    pass_fragments = []
    while len(pass_fragments) < pass_threshold:
        pass_fragments = []
        iteration_counter += 1

        min_error = None
        filter_with_min_error = None
        fragment_with_min_error = None

        for fragment in fragments_to_relax:
            molecule_pass = True
            molecule_error = 0

            mol_filter_min_error = None
            mol_filter_with_min_error = None

            for filter_ in filters:
                distances = fragment["distances"]
                values = fragment[filter_["prop"]]
                is_ok, errors = filter_module.eval_filter_from_definition(
                    filter_, distances, values)

                normalized_error = [filter_["normalize"][index](error)
                                    for index, error in enumerate(errors)]

                molecule_pass &= is_ok
                filter_error = sum([e * e for e in normalized_error])
                molecule_error += filter_error
                if not is_ok and (mol_filter_with_min_error is None or
                                  filter_error < mol_filter_min_error):
                    mol_filter_min_error = filter_error
                    mol_filter_with_min_error = filter_

            if molecule_pass:
                pass_fragments.append(fragment)
                continue

            if (min_error is None or min_error > molecule_error) and \
                    mol_filter_with_min_error is not None:
                min_error = molecule_error
                filter_with_min_error = mol_filter_with_min_error
                fragment_with_min_error = fragment

        # It is possible that we have nothing to improve, ie. all molecules
        # pass the filters.
        if fragment_with_min_error is None:
            assert len(fragments_to_relax) == len(pass_fragments)
            break

        output += " iteration: {:>2} ".format(iteration_counter)
        output += " pass: {:>2}  ".format(len(pass_fragments))
        output += " relaxing: {:<38}".format(fragment_with_min_error["name"])
        output += " with error: {:>6.4f}\n".format(min_error)

        output += "    {:<8}  {:<15}  {:<40}  {}\n".format(filter_with_min_error["@type"],
                  filter_with_min_error["prop"],
                  ", ".join(["[{:>4.3f}, {:>4.3f}]".format(*item) for item in filter_with_min_error["intervals"]]),
                  ", ".join(["({} {} {} {})".format(item[0], item[1], item[3], item[4]) for item in filter_with_min_error["keys"]]))

        relax_filter(filter_with_min_error, fragment_with_min_error)
        
        output += "    {:<8}  {:<15}  {:<40}\n\n".format("", "", ", ".join(["[{:>4.3f}, {:>4.3f}]".format(*item)
                 for item in filter_with_min_error["intervals"]]))

    filters_to_remove = []
    for filter_ in filters:
        if (filter_["prop"] == "distances") and \
           (filter_["intervals"][0][1] > 5.0):
                logger.debug("Distance over threshold: %s", filter_)
                filters_to_remove.extend(filter_["keys"])
                continue
        if (filter_["prop"] == "accessibility") and \
           (filter_["intervals"][0][0] < 0.5):
                logger.debug("Accessibility under threshold: %s", filter_)
                filters_to_remove.extend(filter_["keys"])
    if len(filters_to_remove) > 0:
        return orig_filters, None, filters_to_remove, ""

    # Remove normalize functions.
    for filter_ in filters:
        del filter_["normalize"]

    output += "Filters:\n"
    output += " \n  ".join(map(str, filters))
    output += "\n"
    output += "Failed fragments:\n"
    print(output)
    print_details(fragments, filters)

    print("\n" * 4)

    metadata = {
        "iterations": iteration_counter,
        "input_size": len(fragments),
        "threshold": pass_threshold,
        "final_pass": len(pass_fragments),
        "bad_molecules": len(bad_fragments),
        "used_for_relaxation": len(fragments)
    }
    return filters, metadata, None, output


def _remove_filter(filters, to_discard):
    logger.info("Removing nonsense interactions: %s", to_discard)
    new_filters = []
    for fil in filters:
        for td in to_discard:
            if fil["keys"][0] == td:
                break
        else:
            new_filters.append(fil)
    return new_filters


# endregion


def print_details(fragments, filters):
    for fragment in fragments:
        values_dict = {}
        for filter_ in filters:
            if filter_["prop"] in values_dict:
                continue
            values_dict[filter_["prop"]] = fragment[filter_["prop"]]
        filter_module.eval_debug_filters_from_definitions(
            filters, values_dict, fragment["name"])


def _check_filters_dist_access(filters):
    access_keys = set()
    filters_to_use = []
    for fil in filters:
        if fil["prop"] == "accessibility":
            for _key in fil["keys"]:
                access_keys.add(tuple(_key))
            filters_to_use.append(fil)
    for fil in filters:
        if fil["prop"] == "distances":
            if tuple(fil["keys"][0]) in access_keys:
                filters_to_use.append(fil)
    return filters_to_use


if __name__ == "__main__":
    main()
