#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import typing
import collections

from libs import paths
from libs import file_formats
from libs.data_model import McsaEntry, MappingEntry
from libs.configuration import CHEM, TASK


class ActiveSite(object):
    class Residue(object):
        def __init__(self, res_seq, res_name, type_, atoms):
            self.res_seq = res_seq
            self.res_name = res_name
            self.type = type_
            self.atoms = atoms

    def __init__(self):
        self._residues = []

    def residues(self) -> typing.List[Residue]:
        return self._residues

    def residue_by_seq(self, res_seq):
        for residue in self._residues:
            if residue.res_seq == res_seq:
                return residue
        return None


def activity_threshold(res_name, type_, prop_to_asses=None):
    key = (res_name, type_)
    res_function = CHEM["residues"].get(key, None)
    if res_function is None:
        return None
    if prop_to_asses == "accessibility":
        return res_function["accessibility_threshold"]
    return res_function["threshold"]


def activity_tolerance(res_name, type_):
    key = (res_name, type_)
    res_function = CHEM["residues"].get(key, None)
    if res_function is None:
        return None
    return res_function["tolerance"]


def active_site_from_mcsa(mcsa: McsaEntry):
    active_site = ActiveSite()
    for residue in mcsa.residues().first().value():
        res_function, type_ = _function_for_residue(residue)
        if res_function is None:
            message = "Missing interaction description for "
            message += residue["code"] + " " + str(residue["position"]) + " "
            message += str(residue["types"])
            raise RuntimeError(message)
        active_site.residues().append(ActiveSite.Residue(
            residue["position"],
            residue["code"],
            type_,
            res_function["atoms"]
        ))
    return active_site


def _function_for_residue(residue):
    for type_ in residue["types"]:
        key = (residue["code"], type_)
        res_function = CHEM["residues"].get(key, None)
        if res_function is not None:
            return res_function, type_
    return None, None


# TODO Move to file where used -> collect_Data
def all_residues_atoms(pdbqt: file_formats.PdbqtObject,
                       mapping: MappingEntry,
                       active_site: ActiveSite):
    chain = mapping.source_chain().first().value()
    offset = mapping.offset().first().value()
    atom_pos_mapping = _reverse_dict(mapping.fasta_to_mcsa().first().value())

    res_seq_to_residue = {}
    for residue in active_site.residues():
        res_seq = atom_pos_mapping[residue.res_seq] + offset
        res_seq_to_residue[res_seq] = residue

    output = {}
    for atom in pdbqt.atoms():
        if not atom.chain == chain:
            continue
        if atom.res_seq not in res_seq_to_residue:
            continue
        residue = res_seq_to_residue[atom.res_seq]
        if residue not in output:
            output[residue] = [atom]
        else:
            output[residue].append(atom)
    return output


# TODO Move to file where used -> docking_evaluate, ??
def reacting_atoms(pdbqt: file_formats.PdbqtObject,
                   mapping: MappingEntry,
                   active_site: ActiveSite):
    chain = mapping.source_chain().first().value()
    key_to_name = _create_atom_key_to_name(mapping, active_site)

    output = collections.defaultdict(lambda: [])
    used_keys = set()
    for atom in pdbqt.atoms():
        if not atom.chain == chain:
            continue
        key = _create_atom_key(atom)
        name = key_to_name.get(key, None)
        if name is None:
            continue
        used_keys.add(key)
        output[name].append(atom)

    if len(key_to_name.keys() - used_keys) == 0:
        return output

    # See issue #4 for example of failing data.
    error_report = "Chain: {}\n".format(chain)
    error_report += "Filter:\n  "
    error_report += "\n  ".join(map(str, key_to_name))
    error_report += "\nAtoms founds:\n  "
    error_report += "\n  ".join(["{} : {}".format(key, str(len(value)))
                                 for key, value in output.items()])
    error_report += "\n"
    raise RuntimeError("Invalid number of atoms found. \n" + error_report)


def _create_atom_key_to_name(mapping, active_site):
    offset = mapping.offset().first().value()
    atom_pos_mapping = _reverse_dict(mapping.fasta_to_mcsa().first().value())

    key_to_name = {}
    for residue in active_site.residues():
        name = (residue.res_seq, residue.res_name, residue.type)
        for atom in residue.atoms:
            res_seq = atom_pos_mapping[residue.res_seq] + offset
            atom_filter = (res_seq, residue.res_name, atom)
            key_to_name[atom_filter] = name
    return key_to_name


def _create_atom_key(atom):
    # TODO Issue #3
    name = atom.name
    if len(name) == 4:
        """
        Cases where there is reference for HD21 and HD22  but
        the pdbqt file use 1HD2 and 2HD2 format instead.
        Now, we follow the standard PDB residue naming convention 
        for the reference: ftp://ftp.wwpdb.org/pub/pdb/data/monomers
        """
        name = name[1:] + name[0]
    return atom.res_seq, atom.res_name, name


def _reverse_dict(dictionary):
    return {val: int(key) for key, val in dictionary.items()}


def load_reactive_atoms(mapping: MappingEntry, mcsa: McsaEntry):
    protein_path = paths.protein_pdbqt_path(mapping, mcsa)
    protein_pdbqt = file_formats.load_pdbqt_file(protein_path)
    active_site = active_site_from_mcsa(mcsa)
    atoms = reacting_atoms(protein_pdbqt, mapping, active_site)
    return atoms


def compute_docking_box(mapping: MappingEntry, mcsa: McsaEntry):
    protein_path = os.path.join(
        paths.proteins(), paths.file_name_for_mapping(mapping, mcsa) + ".pdbqt")
    pdbqt_file = file_formats.load_pdbqt_file(protein_path)
    active_site = active_site_from_mcsa(mcsa)
    atoms_by_residues = all_residues_atoms(pdbqt_file, mapping, active_site)
    # TODO Check that center is computed correctly.
    center = _atoms_center(atoms_by_residues)
    box = TASK["docking"]["box"]
    return center, box


def _atoms_center(atoms_by_residues):
    center_x, center_y, center_z = [], [], []
    for atoms in atoms_by_residues.values():
        x, y, z = 0, 0, 0
        for item in atoms:
            x += item.x
            y += item.y
            z += item.z
        size = len(atoms)
        center_x.append(x / float(size))
        center_y.append(y / float(size))
        center_z.append(z / float(size))

    def avg(values):
        return sum(values) / float(len(values))

    return avg(center_x), avg(center_y), avg(center_z)
