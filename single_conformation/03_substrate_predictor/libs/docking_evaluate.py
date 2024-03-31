#!/usr/bin/env python
# -*- coding: utf-8 -*-

import typing
import math
import logging
import numpy as np

from libs import file_formats
from libs.data_model import ReadOnlyDatabase
from libs import vector_math as vmath
from libs.configuration import CHEM


class FragmentContainer(object):
    def __init__(self, fragment, molecule, atoms, atoms_name):
        self.fragment = fragment
        self.molecule = molecule
        self.atoms = atoms
        self.atoms_name = atoms_name
        self.filter_atom_dist = None
        self.filter_charge_dist = None


logger = logging.getLogger(__name__)

# region loading docking results


def load_docking_result(pdbqt_path, sdf_path, docking_path):
    docking_result = file_formats.load_pdbqt_models_file(docking_path)

    substrate_sdf = file_formats.load_sdf_file(sdf_path)
    substrate_pdbqt = file_formats.load_pdbqt_file(pdbqt_path)
    file_formats.add_sdf_bonds_to_pdbqt(substrate_sdf, substrate_pdbqt)

    # We assume that the order of atoms is the same.
    for mol in docking_result:
        mol._bonds = substrate_pdbqt._bonds
    return docking_result


def compute_docking_cost():
    from scipy.optimize import curve_fit
    from libs import paths
    import json
    import os

    def _fit_function(x, a, b, c):
        return a * np.exp(b * x) + c

    cost_function = None
    max_torsdofs_times = {}
    for filename in os.listdir(paths.complexes()):
        if filename.endswith("_time.json"):
            json_path = os.path.join(paths.complexes(), filename)
            with open(json_path) as in_stream:
                substrate_description = json.load(in_stream)
            if substrate_description["torsdof"] not in max_torsdofs_times.keys():
                max_torsdofs_times[substrate_description["torsdof"]] = []

            max_torsdofs_times[substrate_description["torsdof"]].append(substrate_description["docking_time"])

    if max_torsdofs_times:
        torsdofs = [float(key) for key in sorted(max_torsdofs_times.keys())]
        max_times = [max(max_torsdofs_times[key]) for key in sorted(max_torsdofs_times.keys())]
        popt, pcov = curve_fit(_fit_function, torsdofs, max_times)

        cost_function = {}
        for i in range(0, 21):
            cost_function[i] = int(round(_fit_function(i, *popt)))

    return cost_function


# endregion

# region fragment extraction


def extract_fragments(fragment_definition, molecule):
    """
    Extract fragments (atoms) that satisfy the definition. Return
    list of lists of atoms that create the fragment. The order of atoms
    is same as in fragment_definition.atoms.

    TODO:
    Correctness of this implementation should be checked. An alternative
    is to use some chemical library with fragment matching function, that
    would however introduce a dependency.
    """
    # print(fragment_definition)
    # print(molecule)
    candidates = collect_candidates(fragment_definition, molecule)
    # print(candidates)
    # print(candidates[0][0].x, candidates[1][0].x)
    candidates = filter_by_bond(fragment_definition, molecule, candidates)
    # print(candidates)
    # quit()
    # print(str(candidates[1][0]), str(candidates[1][1]))
    # print(str(candidates[2][0]))
    # quit()
    min_candidates_count = min(map(len, candidates))
    if min_candidates_count == 0:
        return []
    return collect_fragments(fragment_definition, molecule, candidates)


def collect_candidates(fragment_definition, molecule):
    candidates = [[] for _ in fragment_definition["types"]]
    for atom in molecule.atoms():
        for index in range(len(candidates)):
            candidate_type = fragment_definition["types"][index]
            if are_atom_types_matching(candidate_type, atom.atom_type):
                candidates[index].append(atom)
    return candidates


def are_atom_types_matching(atom_template, atom_type):
    # TODO Add mapping to the knowledge base.
    if atom_template == atom_type:
        return True
    if atom_template == "X":
        return atom_type in {"F", "Cl", "Br", "I", "AT"}
    if atom_type == "OA":
        return atom_template in {"O"}
    if atom_type == "NA":
        return atom_template in {"N"}
    if atom_type == "A":
        return atom_template in {"C"}
    return False


def filter_by_bond(fragment, molecule, candidates):
    output = [[] for _ in candidates]
    for bond in fragment["bonds"]:
        for origin in candidates[bond[0]]:
            for target in candidates[bond[1]]:
                if molecule.contains_bond(origin, target, bond[2]):
                    output[bond[0]].append(origin)
                    output[bond[1]].append(target)
    return output


def collect_fragments(fragment, molecule, candidates):
    """
    We try all combinations and then filter such that does not contains
    required bonds.
    """
    output = []
    for seed in candidates[0]:
        output.extend(select_subsets(candidates, [seed]))
    filtered_output = [item for item in output
                       if do_satisfy_all_bonds(fragment, molecule, item)]
    return filtered_output


def select_subsets(candidates, accumulator):
    if len(accumulator) == len(candidates):
        # We have picked one molecule from each candidate.
        return [accumulator]
    next_index = len(accumulator)
    output = []
    for next_mol in candidates[next_index]:
        output.extend(select_subsets(candidates, accumulator + [next_mol]))
    return output


def do_satisfy_all_bonds(fragment, molecule, atoms):
    for bond in fragment["bonds"]:
        if not molecule.contains_bond(atoms[bond[0]], atoms[bond[1]], bond[2]):
            return False
    return True


# endregion

# region basic fragment evaluation (cheap operations)

def evaluate_fragments_basic(database: ReadOnlyDatabase,
                             fragments: typing.List[FragmentContainer],
                             conformer: int=None):
    for fragment in fragments:
        evaluate_fragment_basic(database, fragment, conformer)


def evaluate_fragment_basic(database: ReadOnlyDatabase,
                            container: FragmentContainer,
                            conformer):
    fragment = container.fragment
    mapping_ref = fragment.mapping
    if conformer is None:
        protein_atoms = load_protein(database, mapping_ref)
    else:
        protein_atoms = load_protein(database, mapping_ref, conformer)
    active_atoms = protein_atoms["reactive"]
    try:
        container.filter_atom_dist = evaluate_atom_distance(
            active_atoms, container)
    except IndexError as ex:
        logger.error("Can't evaluate distance for '%s'", fragment.name)
        raise ex
    container.filter_charge_dist = evaluate_charge_distances(
        active_atoms, container)
    fragment.coulomb_score = coulombs_score(active_atoms, container)
    # TODO Temporary.
    # fragment.props["angle"] = temporary_compute_angle(active_atoms, container)


def load_protein(database, mapping_ref):
    raise NotImplementedError("Implementation must be provided by API user.")


def evaluate_atom_distance(active_site, container: FragmentContainer):
    """
    Each atom in a fragment must be closer then ACTIVITY_CUT_OFF,
    to any active site molecule.
    """
    for atom in container.atoms:
        dist = minimum_atoms_distance(active_site, [atom])
        if dist > CHEM["activity_threshold"]:
            return False
    return True


def evaluate_charge_distances(active_site, container: FragmentContainer):
    neg_fragment = [atom for atom in container.atoms if atom.charge() < 0]
    pos_active_site = [atom for atom in active_site if atom.charge() > 0]
    for neg_atom in neg_fragment:
        dist = minimum_atoms_distance(pos_active_site, [neg_atom])
        if dist > CHEM["activity_threshold"]:
            return False

    pos_fragment = [atom for atom in container.atoms if atom.charge() > 0]
    neg_active_site = [atom for atom in active_site if atom.charge() < 0]
    for pos_atom in pos_fragment:
        dist = minimum_atoms_distance(neg_active_site, [pos_atom])
        if dist > CHEM["activity_threshold"]:
            return False

    return True


def minimum_atoms_distance(left_atoms, right_atoms):
    output = atoms_distance(left_atoms[0], right_atoms[0])
    for left in left_atoms:
        for right in right_atoms:
            output = min(output, atoms_distance(left, right))
    return output


def atoms_distance(left_atom, right_atom):
    distances = [left_atom.x - right_atom.x,
                 left_atom.y - right_atom.y,
                 left_atom.z - right_atom.z]
    return math.sqrt(sum([x * x for x in distances]))


def coulombs_score(active_site, container: FragmentContainer):
    """
    Coulomb's law without the Coulomb's constant.
    """
    total_energy = 0
    for protein_atom in active_site:
        for fragment_atom in container.atoms:
            r = atoms_distance(protein_atom, fragment_atom)
            energy = (protein_atom.charge() * fragment_atom.charge()) / (r * r)
            total_energy += energy
    return total_energy


# endregion

# noinspection PyPep8Naming
def temporary_compute_angle(active_atoms, fragment):
    C = None
    X = None
    for atom in fragment.atoms:
        if atom.atom_type == "C":
            C = vmath.atom_as_vector(atom)
        else:
            X = vmath.atom_as_vector(atom)
    # Select O from ASP that is closer
    ASP = None
    ASP_C_dist = 0
    for atom in active_atoms:
        if not atom.name.startswith("OD"):
            continue
        atom_vector = vmath.atom_as_vector(atom)
        C_atom_dist = vmath.vector_size(vmath.vector_diff(C, atom_vector))
        if ASP is None:
            ASP = atom_vector
            ASP_C_dist = C_atom_dist
        elif C_atom_dist < ASP_C_dist:
            ASP = atom_vector
            ASP_C_dist = C_atom_dist

    ASP_C = vmath.vector_diff(ASP, C)
    X_C = vmath.vector_diff(X, C)

    cos_alpha = vmath.vector_dot_product(ASP_C, X_C) / (vmath.vector_size(ASP_C) * vmath.vector_size(X_C))
    return math.degrees(math.acos(cos_alpha))


# region fragment evaluation (expensive operations)

def evaluate_fragments_expensive(
        database, containers: typing.List[FragmentContainer], conformer=None):
    for container in containers:
        mapping_ref = container.fragment.mapping
        if conformer is None:
            protein = load_protein(database, mapping_ref)
        else:
            protein = load_protein(database, mapping_ref, conformer)
        distances = compute_distances(container, protein["reactive_by_name"])
        round_values(distances, 3)
        container.fragment.distances = distances


def compute_distances(container: FragmentContainer, atoms_by_name):
    output = []
    for frag_name, frag_atom in zip(container.atoms_name, container.atoms):
        for name, atoms in atoms_by_name.items():
            for atom in atoms:
                name = protein_fragment_interaction_key(
                    name[0], name[1], name[2], atom.name, frag_name)
                output.append({
                    "name": name,
                    "val": atoms_distance(frag_atom, atom)
                })
    return output


def round_values(values_list, ndigits):
    for item in values_list:
        item["val"] = round(item["val"], ndigits=ndigits)


# endregion

def protein_fragment_interaction_key(
        res_seq, res_name, res_type, res_atom_name, frag_atom_name):
    return res_seq, res_name, res_type, res_atom_name, frag_atom_name


def filter_out_bad(containers):
    return [container for container in containers
            if container.fragment.vina_score < 0
            and container.filter_atom_dist
            and container.filter_charge_dist]
