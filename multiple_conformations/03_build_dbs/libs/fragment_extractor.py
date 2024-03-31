#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections

from libs import file_formats


def extract_fragment(sdf_path):
    sdf = file_formats.load_sdf_file(sdf_path)
    fragment = {
        "types": [],
        "names": [],
        "bonds": []
    }
    
    index_by_type = collections.defaultdict(int)
    index_to_atom = {}
    
    index = 0
    for atom in sdf.atoms():
        index += 1
        if atom.atom_type == "R":
            continue
        index_to_atom[index] = len(fragment["types"])
        fragment["types"].append(atom.atom_type)
        index_by_type[atom.atom_type] = index_by_type[atom.atom_type] + 1
        fragment["names"].append(atom.atom_type + str(index_by_type[atom.atom_type]))
    print(index_to_atom, "<- itoa")
    
    for bond in sdf.bonds():
        #print(bond.origin, "<- origin")
        #print(bond.target, "<- target")
        if bond.origin not in index_to_atom or bond.target not in index_to_atom:
            continue
        fragment["bonds"].append([index_to_atom[bond.origin],
                                  index_to_atom[bond.target],
                                  int(bond.type)])
    
    return fragment
