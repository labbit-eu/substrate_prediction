#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import json
import logging
import os
import functools
import collections
import gzip
import multiprocessing as mp

from libs import file_formats
from libs import paths
from libs.data_model import ReadOnlyDatabase
from libs import vector_math as vmath
from libs import utils
from libs.configuration import CHEM, TASK, ENV, get_num_cpu, _SPECIAL_RESIDUES
from libs import interactions
from libs.filter import _load_filters
from external import add_hydrogens_to_ligand

MODE = "substrates"


logger = logging.getLogger(__name__)


def add_accessibility_to_file(input_files, output_files):

    database = ReadOnlyDatabase(paths.database_directory())
    if MODE == "substrates":
        gzipped = False
    else:
        gzipped = True
    for input_file, output_file in zip(input_files, output_files):
        logger.info("Adding accessibility '%s' -> '%s'", input_file, output_file)
        logger.info("Loading fragments ...")
        fragments = utils.load_json_lines(input_file, gzipped)
        logger.info("Fragments count: %d", len(fragments))
        if len(fragments) == 0:
            logger.info("Skipping file.", len(fragments))
            return
        logger.info("Collecting molecules for protonation ...")
        # hydrogen_dir = os.path.join(paths.temp(), _get_temp_directory())
        # prepare_sdf_for_protonation(database, fragments, hydrogen_dir)
        prepare_sdf_for_protonation(database, fragments)
        logger.info("Executing PyMol script ...")
        add_hydrogens_to_ligand.add_hydrogens_batch(paths.temp(), MODE)

        if MODE == "substrates":
            candidate_filters_names = _create_candidate_filters_names_substrates(database, fragments)
        else:
            candidate_filters_names = _extract_filters_names_candidates(input_file)

        logger.info("Computing accessibility ...")
        pool = mp.Pool(processes=get_num_cpu(ENV["num_parallel_cpu"]))
        processing = []
        new_fragments = []

        for fragments_set in utils.split_fragments_to_sets(fragments):
            processing.append(pool.apply_async(_add_accessibility_to_set_of_fragments,
                                               args=(database, fragments_set, candidate_filters_names)))
        for p in processing:
            new_fragments += p.get()

        pool.close()
        logger.info("In total, %d candidates processes.", len(new_fragments))

        _save_fragments(new_fragments, output_file)


def _add_accessibility_to_set_of_fragments(database, fragments_set, candidate_names):

    for fragment in fragments_set:
        _conf = fragment["ref"]["mapping"].split(",")[1]
        _conf = _conf.split("-")[0]
        _conf = _conf.split("_")[-1][1:]
        conformer = int(_conf)
        _add_accessibility_to_fragment(database, fragment, candidate_names, conformer)

    logger.info("Computed accessibility for %d fragments", len(fragments_set))

    return fragments_set


def _add_accessibility_to_fragment(database, fragment, candidate_names, conformer):
    protein = load_protein_atoms(database, fragment["ref"]["mapping"], conformer)

    interactions = _interacting_pairs_from_names(fragment, protein, candidate_names)

    # Load molecule with hydrogens and replace molecules in reactive molecules.
    molecule_with_hydrogens = _load_molecule_with_hydrogens(fragment)
    _replace_candidate_molecules(interactions, molecule_with_hydrogens, fragment)

    visualisation = []

    accessibility = []
    for interaction in interactions:
        name = (*interaction["site_name"], interaction["site_atom"].name,
                interaction["frag_name"])

        val, source, targets, outcomes = _atom_to_atom_accessibility(
            interaction["frag_atom"], interaction["site_atom"],
            molecule_with_hydrogens.atoms())

        _update_visualisation(visualisation, source, targets, outcomes,
                              interaction["site_name"][1])

        if (interaction["site_name"][1] in _SPECIAL_RESIDUES) and (interaction["site_atom"].name != "NZ"):
            name = (*interaction["site_name"], _SPECIAL_RESIDUES[name[1]][interaction["site_atom"].name], interaction["frag_name"])
        accessibility.append({
            "name": name,
            "val": val
        })

    if TASK[MODE]["accessibility"]["save_visualisation"]:
        _save_visualisation(visualisation, fragment)

    fragment["accessibility"] = accessibility


def _get_temp_directory(conformer=None):
    if MODE == "substrates":
        return os.path.join("sdf-hydrogens", "substrates")
    else:
        if conformer is None:
            return "sdf-hydrogens"
        else:
            return os.path.join("conformer_c{}".format(conformer), "sdf-hydrogens")


def prepare_sdf_for_protonation(database, fragments, directory):
    raise NotImplementedError("Provide implementation of this method.")


def _save_fragments(fragments, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if MODE == "substrates":
        with open(output_path, "w") as out_stream:
            for fragment in fragments:
                json.dump([fragment], out_stream)
                out_stream.write("\n")
    else:
        with gzip.open(output_path, "wt") as out_stream:
            for fragment in fragments:
                json.dump([fragment], out_stream)
                out_stream.write("\n")


def load_protein_atoms(database, mapping_ref):
    raise NotImplementedError("Provide implementation of this method.")


def load_molecule(name):
    raise NotImplementedError("Provide implementation of this method.")


# region candidates for accessibility (interaction)
def _extract_filters_names_candidates(input_path):
    # as here we apply filters, we should just add all distances in particular filter as relevant
    name = os.path.basename(input_path)
    name = name[:name.find(".json")]
    filter_path = os.path.join(paths.filters(), name + ".json")
    logger.info("Extracting relevant distances for computation from the respective filter %s", filter_path)
    filters, used_props = _load_filters(filter_path, {})
    candidate_names = []
    for filter_ in filters:
        filter_name = filter_["keys"]
        if filter_name[0] not in candidate_names:
            candidate_names += filter_name
    return candidate_names


def _create_candidate_filters_names_substrates(database, fragments):
    logger.info("Selecting candidates for computation ...")
    distances = _select_relevant_distances(database, fragments)  # those that have opposite charges on interacting atoms
    relevant_distances = interactions.filter_distances_by_type(
        distances)  # those that have median < threshold in CHEM["residues"]

    return [*relevant_distances.keys()]


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


def _interacting_pairs_from_names(fragment, protein, names):
    candidates = []
    reactive_atoms = protein["reactive_by_name"]
    for frag_atom in fragment["atoms"]:
        for name, atoms in reactive_atoms.items():
            for atom in atoms:
                if name[1] in _SPECIAL_RESIDUES and (atom.name != "NZ"):
                    for dist in fragment["distances"]:
                        if dist["name"][3] == "NZ": continue
                        _dname = dist["name"]
                        if name[0] == _dname[0] and name[1] == _dname[1] and \
                           name[2] == _dname[2] and frag_atom["name"] == _dname[4]:
                            if dist["old_name"]==atom.name:
                                candidate_name = (*name, _SPECIAL_RESIDUES[name[1]][atom.name],
                                                  frag_atom["name"])
                                if candidate_name in names:
                                    candidates.append({
                                        "frag_atom": frag_atom,
                                        "site_name": name,
                                        "site_atom": atom
                                    })
                                break
                else:
                    candidate_name = (*name, atom.name, frag_atom["name"])
                    if candidate_name in names:
                        candidates.append({
                            "frag_atom": frag_atom,
                            "site_name": name,
                            "site_atom": atom
                        })
    return candidates


# endregion

# region adding hydrogens

def _load_molecule_with_hydrogens(fragment):

    name = "{}-{}".format(fragment["name"], fragment["model_index"])
    if MODE == "substrates":
        hydrogen_dir = os.path.join(paths.temp(), _get_temp_directory())
        output_sdf_path = os.path.join(hydrogen_dir, name + "_h.sdf")
    else:
        _conf = fragment["ref"]["mapping"].split(",")[1]
        _conf = _conf.split("-")[0]
        _conf = _conf.split("_")[-1][1:]
        conformer = int(_conf)
        hydrogen_dir = os.path.join(paths.temp(), _get_temp_directory(conformer))
        output_sdf_path = os.path.join(hydrogen_dir, fragment["subdir_name"], name + "_h.sdf")

    if not os.path.exists(output_sdf_path):
        raise RuntimeError("Missing molecule with hydrogens for: " + name)

    return _load_sdf_with_hydrogens(output_sdf_path)


@functools.lru_cache(maxsize=64)
def _load_sdf_with_hydrogens(sdf_path):
    return file_formats.load_sdf_file(sdf_path)


def _replace_candidate_molecules(candidates, molecule_with_hydrogens, fragment):
    """
    Load molecule with hydrogens and replace fragment atoms in candidates
    with corresponding in molecule with hydrogens.
    """
    for candidate in candidates:
        atom = candidate["frag_atom"]
        candidate["frag_name"] = candidate["frag_atom"]["name"]
        candidate["frag_atom"] = _select_atom_by_position(
            [atom["x"], atom["y"], atom["z"]],
            molecule_with_hydrogens.atoms())
        if candidate["frag_atom"] is None:
            print("ATOM")
            print(atom)
            print("MOLECULES")
            print("\n".join(map(str, molecule_with_hydrogens.atoms())))
            print(fragment)
            raise RuntimeError("Can't find atom by position.")


def _select_atom_by_position(pos, atoms):
    for atom in atoms:
        if atom.x == pos[0] and atom.y == pos[1] and atom.z == pos[2]:
            return atom
    return None


# endregion

# region accessibility


def _atom_to_atom_accessibility(fragment_atom, filter_atom, collision_atoms):
    fragment_radius = CHEM["atom_radius"][fragment_atom.atom_type]

    candidates = _select_candidates_for_collision(
        filter_atom, fragment_atom, fragment_radius, collision_atoms)

    fragment_vec = vmath.atom_as_vector(fragment_atom)
    filter_vec = vmath.atom_as_vector(filter_atom)
    source_radius = CHEM["atom_radius"][filter_atom.atom_type]

    fragment_to_filter = vmath.vector_normalized(
        vmath.vector_diff(filter_vec, fragment_vec))

    # Take the closes point on the filter atom.
    # source_point = vmath.vector_add(
    #     vmath.multiply_by_scalar(filter_to_fragment, filter_radius),
    #     filter_vec)
    source_point = filter_vec

    target_points, target_weights = _sample_sphere(
        fragment_to_filter, fragment_radius, fragment_vec)

    outcomes = []
    accessible_surface_est = 0.0
    for target_point, target_weights in zip(target_points, target_weights):
        # If the atoms are too close, the target_point may
        # be "blocked" by the target atom it self. Ie. the target point
        # is not directly visible from source.

        # TODO This return None sometimes, that should not happen.
        # Possible cause can be mathematical instability.
        _, distance_to_target_surface = vmath.do_sphere_obstruct_connection(
            source_point, target_point, fragment_vec, fragment_radius)

        for atom in candidates:
            point = vmath.atom_as_vector(atom)
            radius = CHEM["atom_radius"][atom.atom_type]
            collide, distance = vmath.do_sphere_obstruct_connection(
                source_point, target_point, point, radius,
                ignore_in_radius=source_radius)
            if collide:
                if distance_to_target_surface is not None and \
                        distance_to_target_surface < distance:
                    continue
                outcomes.append(False)
                break
        else:
            accessible_surface_est += target_weights
            outcomes.append(True)

    accessible_surface_est = round(accessible_surface_est, 3)

    return accessible_surface_est, source_point, target_points, outcomes


def _select_candidates_for_collision(
        filter_atom, fragment_atom, fragment_radius, collision_atoms):
    """
    Use distance to line between filter_atom and fragment_atom.
    This also includes atoms that are behind filter or fragment.
    """

    fragment_vec = vmath.atom_as_vector(fragment_atom)
    filter_vec = vmath.atom_as_vector(filter_atom)

    candidates = []
    for atom in collision_atoms:
        if atom == fragment_atom or atom == filter_atom:
            continue

        # Check whether atom is in the proximity of line between
        # fragment and active site atom.
        distance = vmath.point_to_line_distance(
            fragment_vec, filter_vec, vmath.atom_as_vector(atom))
        atom_radius = CHEM["atom_radius"][atom.atom_type]
        if distance > (fragment_radius + atom_radius):
            continue

        candidates.append(atom)
    return candidates


def _sample_sphere(direction, radius, center):
    """
    Sample sphere (center, radius) around given direction
    (coming from center of sphere).

    This is done by rotating the template ~ UNIT_SPHERE_POSITIONS_WITH_WEIGHT.
    """
    sphere_template = _create_sphere_template()

    direction = vmath.vector_normalized(direction)

    weights = []
    sphere_direction = sphere_template[0][0]
    rotation_axis = vmath.vector_product(direction, sphere_direction)
    rotation_axis = vmath.vector_normalized(rotation_axis)

    angle = math.acos(vmath.vector_dot_product(direction, sphere_direction))
    rotation_matrix = vmath.rotation_matrix_from_vector_and_angle(
        rotation_axis, angle)

    points = []
    for weighted_point in sphere_template:
        vector = vmath.multiply_by_matrix(
            weighted_point[0], rotation_matrix)

        scaled_vector = vmath.multiply_by_scalar(vector, radius)

        points.append(vmath.vector_add(scaled_vector, center))
        weights.append(weighted_point[1])

    return points, weights


@functools.lru_cache(maxsize=1)
def _create_sphere_template():
    def float_range(start, end, step):
        output_values = []
        val = start
        while val < end:
            output_values.append(val)
            val += step
        return output_values

    output = []
    for level in CHEM["accessibility"]["sphere_definition"]:
        alpha = level[0]
        beta_interval = float_range(*level[1])
        point_weight = level[2] / len(beta_interval)
        for beta in beta_interval:
            point = [math.cos(beta) * math.sin(alpha),
                     math.cos(alpha),
                     math.sin(beta) * math.sin(alpha)]
            output.append([[round(x, 3) for x in point], point_weight])
    return output


# endregion

# region visualisation


def _update_visualisation(visualisation, source, targets, outcomes, res_name):
    if len(visualisation) == 0:
        res_seq = 0
    else:
        res_seq = visualisation[-1][0] + 1

    visualisation.append((2, source, res_seq, res_name))
    for point, outcome in zip(targets, outcomes):
        visualisation.append((int(outcome), point))


def _save_visualisation(visualisation, fragment):
    path = paths.accessibility_visualisation(fragment)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".json", "w") as out_stream:
        json.dump({
            "visualisation": visualisation,
            "fragment": [atom["serial"] for atom in fragment["atoms"]]
        }, out_stream)

# endregion
