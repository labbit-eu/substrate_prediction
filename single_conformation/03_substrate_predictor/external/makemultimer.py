#!/usr/bin/python
"""
This program expands molecular coordinates in pdb files using the BIOMT
transformation matrix. Annotations other than coordinates are currently
discarded. Support for secondary structure annotation may or may not be added at
a later time. Most viewers do a pretty good job at inferring secondary
structure anyway, however.

Send comments or bug reports to mpalmer at uwaterloo dot ca.
---
Copyright Michael Palmer, University of Waterloo.
Free for use in any way you see fit but you have to leave this notice
in place and in effect in derived works.

This is Python3 compatible modified version that does not rename the chains in monomers.

"""

import sys
import string

# this phrase precedes a matrix
group_start_marker = "APPLY THE FOLLOWING TO CHAINS:"
chain_declaration_marker = "CHAINS:"  # applies to first and subsequent lines
url = "watcut.uwaterloo.ca/makemultimer"


def table_format(title_list, data_lists, row_joiner="  "):
    """
    generic table formatting, adapted to this simple case.
    first row contains titles, all others contain values.
    """
    # first, make sure all numeric entries get properly converted
    raw_rows = data_lists
    row_lists = [title_list[:]]

    for row in raw_rows:
        row_lists.append([str(field) for field in row])

    # bring all lists to equal length
    ml = max([len(row) for row in row_lists])
    for row in row_lists:
        if len(row) < ml:
            row.extend([""] * (ml - len(row)))

    widths = []  # determine required width for each column

    for field in range(ml):
        width = max([len(row[field]) for row in row_lists]) + len(row_joiner)
        widths.append(width)

    joined_rows = []

    for row in row_lists:
        nr = [(s + row_joiner).rjust(widths[i]) for i, s in enumerate(row)]
        joined_rows.append("".join(nr).rstrip())  # remove trailing spaces

    maxl = max([len(row) for row in joined_rows])
    s = "-" * maxl
    joined_rows.append(s)
    joined_rows.insert(1, s)
    joined_rows.insert(0, s)

    return joined_rows


class PdbError(Exception):
    pass


class Atom(object):
    """
    ATOM    687  CA  SER B  92      26.506  23.631  11.670  1.00 63.47           C
    HETATM 6747  K     K G   1      34.346  68.181   0.000  0.50133.54           K
    prepare a template from this for output with transformed coordinates.
    """
    def __init__(self, rawline):
        self.rawline = rawline
        self.orig_number = self.number = int(rawline[6:11])
        self.chain = rawline[21]
        self.residue = int(rawline[22:26])

        coords = []
        for i in range(30, 54, 8):
            coords.append(float(rawline[i:i+8]))

        self.coords = coords

        info0 = rawline[:6]
        info1 = rawline[11:21]
        info2 = rawline[26:30]
        info3 = rawline[54:]

        self.descriptor_template = "%s%%5d%s%%s%%4d%s" % (info0, info1, info2)
        self.coordinate_template = "%%8.3f%%8.3f%%8.3f%s" % info3

    def transformed(self, *coords):
        """
        fill in new coordinates, leave descriptor template blank -
        will be filled in later.
        """
        return self.descriptor_template + self.coordinate_template % coords

    def __str__(self):
        return self.rawline

    __repr__ = __str__


class ReplicationGroup(object):
    """
    one biomolecule may contain several matrices that apply to
    separate sets of chains. This class represents one matrix
    and the chains it is to be applied to.
    We do the coordinate replication here, but leave the
    naming and numbering to the BioMolecule class.
    """
    def __init__(self, bm_lines, original_chains, options):
        self.original_chains = original_chains  # reference to all chains in the pdb
        self.options = options
        self.source_chains = set()
        self.replicated_chains = {}
        self.matrices = []

        matrix_lines = []

        for line in bm_lines:
            found = line.find(chain_declaration_marker)
            if found > -1:
                raw_chains = line[(found + len(chain_declaration_marker)):]
                sc = [r.strip(",") for r in raw_chains.split()]
                self.source_chains.update(sc)

            elif line.startswith("BIOMT"):
                matrix_lines.append(line)

        # now, process matrix lines
        while matrix_lines:
            first, matrix_lines = matrix_lines[:3], matrix_lines[3:]
            new_matrix = []
            for j in range(3):
                ml = first[j][11:]
                frags = ml.split()
                new_matrix.append([float(x) for x in ml.split()])
            self.matrices.append(new_matrix)

        for chain in sorted(list(self.source_chains)):
            self.replicated_chains[chain] = self._replicate_chain(chain)

    def _replicate_chain(self, chain):
        """
        apply all applicable transformations to one chain and return all
        resulting copies
        """
        atoms = self.original_chains[chain]
        replicated = []

        # for i, matrix in enumerate(self.matrices):
        for matrix in self.matrices:
            transformed_atoms = []

            for atom in atoms:
                coords = []
                x, y, z = atom.coords
                for a, b, c, d in matrix:
                    coords.append(a * x + b * y + c * z + d)

                t = atom.transformed(*coords)
                transformed_atoms.append((t, atom.residue))

            replicated.append(transformed_atoms)

        return replicated

    def __len__(self):
        """
        calculate the number of chains that will be present after replication.
        """
        return len(self.matrices) * len(self.source_chains)


class BioMolecule(object):
    """
    represents one biomolecule. since the main action is in the ReplicationGroup
    class, the parsing here is rather dumb. the main job here is to delegate
    to one or more ReplicationGroups and later to merge their output as needed.
    """
    title_template = "Multimer expanded from BIOMT matrix in pdb file %s"
    # header_line_template = "chain: %s orig. chain: %s  residues: %4d-%4d  atoms: %5d-%5d"

    chain_titles = ["Chain", "original", "1st resid.",
                    "last resid.", "1st atom", "last atom"]

    def __init__(self, bm_lines, all_chains, options):
        """
        all_chains is a reference to the dictionary with all
        chains available in the entire pdb file, from which
        we will extract what is necessary.
        """
        self.all_chains = all_chains
        self.options = options
        self.replication_groups = []
        self.bm_lines = bm_lines
        self.collected_chains = []
        self.overflow_warnings = set()

        in_group = False
        group = []

        for line in bm_lines:
            if not line.strip():
                continue

            if group_start_marker in line:
                in_group = True
                if group:  # finish up preceding group if present
                    self.replication_groups.append(
                        ReplicationGroup(group, self.all_chains, self.options))
                    group = []

            if in_group:
                group.append(line)

        if group:  # finish up last group
            self.replication_groups.append(
                 ReplicationGroup(group, self.all_chains, self.options))

    def collate(self):
        """
        collate all replicated chains into final output
        apply chain renaming and residue renumbering as
        requested
        """
        # first, assign available letters to chains.
        orig_chains = set()
        for rg in self.replication_groups:
            orig_chains.update(rg.replicated_chains.keys())

        orig_chains = sorted(list(orig_chains))

        # make a generator for each original chain
        # that will supply us with the next chain name
        # and residue offset
        def chain_position_generator(orig_pos):
            chain_list = list(string.ascii_uppercase[orig_pos::len(orig_chains)])
            chain_counter = 0

            while True:
                residue_offset = 0
                for x in range(self.options["renamechains"]):
                    yield chain_list[chain_counter], residue_offset
                    residue_offset += self.options["renumberresidues"]
                chain_counter += 1
                if chain_counter == len(chain_list):
                    self.overflow_warnings.add("Chain name overflow when replicating chain %s" % chain_list[0])
                    chain_counter -= 1  # revert to last available letter

        chain_store = {}
        for pos, chain in enumerate(orig_chains):
            chain_store[chain] = chain_position_generator(pos)
        
        # should we, or shouldn't we renumber atoms uniquely? we will.
        atom_numbers = dict.fromkeys(string.ascii_uppercase, 0)

        for rg in self.replication_groups:
            for old_chain, chains in sorted(rg.replicated_chains.items()):
                for chain in chains:
                    new_chain, residue_offset = next(chain_store[old_chain])
                    # do not rename chains in monomers
                    if len(orig_chains) == 1:
                        new_chain = old_chain

                    lst = []
                    atom_offset = atom_numbers[new_chain]
                    atom_counter = 0

                    for atom, res_no in chain:
                        # calculate atom number
                        atom_counter += 1
                        new_atom_no = atom_offset + atom_counter
                        if new_atom_no == 100000:
                            self.overflow_warnings.add("Atom number overflow when replicating chain %s" % old_chain)
                            new_atom_no = 1

                        # adjust offset of res_no
                        rn_corrected = res_no + residue_offset
                        if rn_corrected > 9999:
                            self.overflow_warnings.add("Residue number overflow when replicating chain %s" % old_chain)
                            rn_corrected = max(1, rn_corrected - 10000)

                        # fill in the atom annotation
                        stuff = (new_atom_no, new_chain, rn_corrected)
                        filled = atom % (new_atom_no, new_chain, rn_corrected)
                        lst.append((new_atom_no, rn_corrected, filled))

                    self.collected_chains.append((old_chain, new_chain, lst[0][1], lst[-1][1], lst[0][0],
                                                  lst[-1][0], [l[2] for l in lst]))

                    # adjust offset for next time around
                    rounded = int(round(len(chain), -3))
                    if rounded < len(chain):
                        rounded += 1000
                    atom_numbers[new_chain] += rounded

    def output(self, filename="stuff"):
        """
        return our collected results in one big string. This will
        be the main part of the output pdb file.
        """
        self.collate()

        header = [self.title_template % filename]
        header.append("by MakeMultimer.py (%s)" % url)
        header.append("")

        if self.overflow_warnings:
            header.append("The following errors occurred with naming and numbering:")
            header.extend(sorted(list(self.overflow_warnings)))
            header.append("")

        atoms_out = []
        chain_data_lines = []

        for old, new, first_res, last_res, first_atom, last_atom, atoms in self.collected_chains:
            chain_data_lines.append((new, old, first_res, last_res, first_atom, last_atom))
            atoms_out.extend(atoms)

        formatted = table_format(self.chain_titles, chain_data_lines)

        header.extend(formatted)

        header.append("")
        header.append("BIOMT instructions applied:")
        header.extend(self.bm_lines)
        header.append("")

        out = ["REMARK  " + h for h in header] + atoms_out

        return "\n".join(out)


class PdbReplicator(object):
    """
    responsible for slicing up one pdb file into biomolecules.
    the actual replication occurs elsewhere.
    """
    # backbone atoms for proteins and DNA/RNA
    backbone_atoms = "C CA O N P O3' O5' C3' C4' C5'".split()  # don't need whole ribose

    def __init__(self, pdb_string, options):
        # generator for next running atom number - used across all chains and hetatms
        self.pdb_lines = pdb_string.splitlines()
        self.options = options
        # if we have 0, it means indeed never rename
        self.options["renamechains"] = self.options["renamechains"] or sys.maxint

        self.originalchains = self.parse_molecule()
        self.output = []

        self.biomolecules = self.parse_biomt()

    def parse_biomt(self):
        """
        carve up the file according to biomolecules, which are defined
        in the REMARK 350 lines.
        """
        bm_lines = [l[10:].strip() for l in self.pdb_lines if l.startswith("REMARK 350")]

        if not bm_lines:    # this file doesn't have any biomt instructions for us.
            raise PdbError("input file does not contain any BIOMT instructions")

        # OK, try to carve up the biomolecules
        biomolecules = []

        in_group = False
        group = []

        bm_marker = "BIOMOLECULE:"

        for line in bm_lines:
            if line.startswith(bm_marker):  # start new group ?
                                            # somehow parsing this again seems redundant.
                                            # need to figure it out some more.
                in_group = True
                if group:
                    biomolecules.append(
                        BioMolecule(group, self.originalchains, self.options))
                    group = []
            if in_group and line.strip():
                group.append(line)

        if group:
            biomolecules.append(BioMolecule(group, self.originalchains, self.options))

        return biomolecules

    def test_atom(self, atom_line):
        """
        determine whether or not a regular atom should be
        included.
        - first 5 characters are "ATOM " and columns 15-17 are "CA ".
        - first 5 characters are "ATOM " and columns 15-17 are " P ".
        these numbers are actually false - it is 14-16
        """
        if not self.options["backbone"]:
            return True
        else:
            atom_type = atom_line[13:16].upper().strip()
            return atom_type in self.backbone_atoms

    def test_hetatm(self, hetatm_line):
        """
        determine whether or not a hetatm line should be kept.
        we could extend this to filter out water for example.
        Two conditions:
        - chain identifier must be present
        - depending on nowater option, water will be stripped
        """
        if self.options["nohetatm"]:
            return False

        orig_chain = hetatm_line[21:22].strip()
        if not orig_chain:
            return False

        hetatm_type = hetatm_line[17:20].upper()
        if self.options["nowater"] and hetatm_type in ("HOH", "DOD", "WAT"):
            return False

        return True

    def parse_molecule(self):
        """
        extract all chains from a pdb file. Problem: HETATM records may or may
        not have a chain identifiers.

        we will look at handling those unnamed hetero atoms again later, once
        we get a better feel for the entire thing.
        """
        atom_lines = []

        for line in self.pdb_lines:
            if line.startswith("ATOM") and self.test_atom(line):
                atom_lines.append(line)
            elif line.startswith("HETATM") and self.test_hetatm(line):
                atom_lines.append(line)

        chains = {}
        for al in atom_lines:
            atom = Atom(al)
            chains.setdefault(atom.chain, []).append(atom)

        return chains


