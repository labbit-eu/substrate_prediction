#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import statistics

from libs import active_site as active_site_lib
from libs.configuration import _SPECIAL_RESIDUES


def select_reactive_distances(fragment, molecule, protein):
    # (108, 'ASP', 'covalent catalysis', 'OD2', 'C1') : distance
    distances = collections.defaultdict(list)
    fragment_atoms = _fragment_reactive_atoms(fragment, molecule)
    protein_atoms = _protein_reactive_atoms(protein)
    for distance in fragment["distances"]:
        res_seq, res_name, res_type, res_atom, frag_atom = distance["name"]
        if (res_name in _SPECIAL_RESIDUES) and (res_atom != "NZ"):
            res_atom = distance["old_name"]
        if not _have_opposite_charge(
                fragment_atoms[frag_atom],
                protein_atoms[(res_seq, res_name, res_type, res_atom)]):
            continue
        distances[tuple(distance["name"])].append(distance["val"])
    return distances


def _fragment_reactive_atoms(fragment, molecule):
    fragment_atoms = {}
    for atom in fragment["atoms"]:
        # TODO Get atoms by serial instead of position?
        mol_atom = _select_atom_by_position(
            [atom["x"], atom["y"], atom["z"]],
            molecule.atoms())
        fragment_atoms[atom["name"]] = mol_atom
    return fragment_atoms


def _protein_reactive_atoms(protein):
    protein_atoms = {}
    for name, atoms in protein["reactive_by_name"].items():
        for atom in atoms:
            protein_atoms[(*name, atom.name)] = atom
    return protein_atoms


def _select_atom_by_position(pos, atoms):
    for atom in atoms:
        if atom.x == pos[0] and atom.y == pos[1] and atom.z == pos[2]:
            return atom
    return None


def _have_opposite_charge(left_atom, right_atom):
    return left_atom.partial_chrg * right_atom.partial_chrg < 0


def add_distances(source, destination):
    for key, value in source.items():
        destination[key].extend(value)


def filter_distances_by_type(distances):
    """
    Used interaction specific threshold.
    """
    output = {}
    for name, values in distances.items():
        res_seq, res_name, res_type, res_atom, frag_atom = name
        threshold = active_site_lib.activity_threshold(res_name, res_type)
        if statistics.median(values) > threshold:
            continue
        output[name] = values
    return output
