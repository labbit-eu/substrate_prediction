#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import typing
import logging

logger = logging.getLogger(__name__)


class _AtomCollection(object):
    def __init__(self):
        self._atoms = []

    def __str__(self):
        return "\n".join([str(a) for a in self._atoms]) + "\n"

    def atoms(self):
        return self._atoms

    def iter_atoms_by_res_seq(self, res_seq):
        for atom in self.atoms:
            if atom.res_seq == res_seq:
                yield atom


class _Writable(object):
    def write_to_stream(self, stream):
        raise NotImplementedError


# region pdbqt

class PdbqtAtom(object):
    def __init__(self, line):
        self.serial = int(line[7:11])
        self.name = line[12:16].strip()
        self.res_name = line[17:20]
        self.chain = line[21]
        try:
            # This may not be here for Vina results.
            self.res_seq = int(line[22:26])
        except:
            self.res_seq = None
        self.x = float(line[31:38])
        self.y = float(line[39:46])
        self.z = float(line[47:54])
        self.occupancy = float(line[55:60])
        self.temp_factor = float(line[61:66])
        self.partial_chrg = float(line[67:76])
        self.atom_type = line[77:79]

    def __str__(self):
        return "{:>5} {:>3} {:>3} {:>4} {:>6} {:>6} {:>6} {:>6}".format(
            self.serial, self.name, self.res_name, self.res_seq,
            self.x, self.y, self.z, self.partial_chrg)

    def charge(self):
        return self.partial_chrg


class PdbqtBond(object):
    def __init__(self, origin, target):
        self.origin = int(origin)
        self.target = int(target)

    def __str__(self):
        return "{:>5} {:>5}".format(self.origin, self.target)


class PdbqtObject(_AtomCollection, _Writable):
    """
    https://www.mdanalysis.org/docs/documentation_pages/coordinates/PDBQT.html

    Beware that bonds use indexes starting from 1.
    """

    def __init__(self):
        _AtomCollection.__init__(self)
        _Writable.__init__(self)
        self._file_content_lines = []
        self._bonds = None
        self._props = {}

    def __str__(self):
        if self._bonds is None:
            return "ATOMS\n " + \
                   "\n ".join(map(str, self._atoms)) + \
                   "\n"
        else:
            return "ATOMS\n " + \
                   "\n ".join(map(str, self._atoms)) + \
                   "\nBONDS\n " + \
                   "\n ".join(map(str, self._bonds)) + \
                   "\n"

    def _add_file_line(self, line):
        self._file_content_lines.append(line)

    def write_to_stream(self, stream):
        stream.write_lines(self._file_content_lines)

    def contains_bond(self, origin, target):
        """
        True if there is bond between given atoms.
        """
        if self._bonds is None:
            raise Exception("No bonds information found.")
        for bond in self._bonds:
            # Consider both directions of bond.
            if bond.origin == origin.serial and bond.target == target.serial:
                return True
            if bond.origin == target.serial and bond.target == origin.serial:
                return True
        return False

    def bonds(self):
        if self._bonds is None:
            raise Exception("No bonds information found.")
        return self._bonds

    def prop(self, name):
        return self._props[name]


class _PdbqtLoader(object):

    def __init__(self, target):
        self._target = target

    def line(self, line):
        self._target._add_file_line(line)
        line = line.rstrip()
        try:
            if line.startswith("ATOM"):
                self._target._atoms.append(PdbqtAtom(line))
            elif line.startswith("REMARK VINA RESULT:"):
                data = line[len("REMARK VINA RESULT:"):]
                tokens = [token for token in data.split(" ") if token]
                self._target._props["vina"] = float(tokens[0])
        except Exception as ex:
            logger.info("Error processing line: %s", line)
            raise ex


def load_pdbqt_file(path):
    pdbqt = PdbqtObject()
    loader = _PdbqtLoader(pdbqt)
    with open(path) as input_stream:
        for line in input_stream:
            loader.line(line)
    return pdbqt


def load_pdbqt_models_file(path: str) -> typing.List[PdbqtObject]:
    output = []
    loader = None
    with open(path) as input_stream:
        for line in input_stream:
            if line.startswith("MODEL"):
                pdbqt = PdbqtObject()
                loader = _PdbqtLoader(pdbqt)
                output.append(pdbqt)
            elif line.startswith("ENDMDL"):
                loader = None
            elif loader is not None:
                loader.line(line)

    return output


def save_to_file(path, writable):
    if not isinstance(writable, _Writable):
        raise ValueError("Object must be writable.")

    with open(path, "w") as output_stream:
        writable.write_to_stream(output_stream)


# endregion

# region pdb

class PdbAtom(object):
    def __init__(self, line):
        self.serial = int(line[7:11])
        self.name = line[12:16].strip()
        self.res_name = line[17:20]
        self.chain = line[21]
        self.res_seq = int(line[22:26])
        self.x = float(line[31:38])
        self.y = float(line[39:46])
        self.z = float(line[47:54])

    def __str__(self):
        return "{:>5} {:>3} {:>3} {:>4} {:>6} {:>6} {:>6}".format(
            self.serial, self.name, self.res_name, self.res_seq,
            self.x, self.y, self.z)


class PdbObject(_AtomCollection):
    """
    http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    """

    def __init__(self):
        _AtomCollection.__init__(self)


def load_pdb_file(path):
    target = PdbObject()
    with open(path) as input_stream:
        for line in input_stream:
            line = line.rstrip()
            if line.startswith("ATOM"):
                target._atoms.append(PdbAtom(line))
    return target


# endregion

# region mol2

class Mol2Atom(object):
    def __init__(self, id_, name, x, y, z, type_):
        self.id = int(id_)
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.type = type_

    def __str__(self):
        return "{:>3} {:>5} {:>6} {:>6} {:>6}".format(
            self.id, self.name,
            self.x, self.y, self.z)


class Mol2Bond(object):
    def __init__(self, id_, origin, target, type_):
        self.id = int(id_)
        self.origin = int(origin)
        self.target = int(target)
        self.type = int(type_)

    def __str__(self):
        return "{:>3} {:>3} {:>3} {:>2}".format(
            self.id, self.origin, self.target, self.type)


class Mol2Molecule(object):
    """
    http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
    """

    def __init__(self):
        self._atoms = []
        self._bonds = []

    def atoms(self):
        return self._atoms

    def bonds(self):
        return self._bonds

    def __str__(self):
        return "ATOMS\n " + \
               "\n ".join(map(str, self._atoms)) + \
               "\nBONDS\n " + \
               "\n ".join(map(str, self._bonds)) + \
               "\n"


def load_mol2(path):
    molecule = Mol2Molecule()
    section = None
    with open(path) as input_stream:
        for line in input_stream:
            line = line.rstrip()
            if line.startswith("@<TRIPOS>"):
                section = line[9:]
            else:
                _load_mol2_section(molecule, section, line)
    return molecule


def _load_mol2_section(molecule, section, line):
    tokens = [x for x in line.rstrip().split(" ") if x]
    if section == "ATOM":
        atom = Mol2Atom(tokens[0], tokens[1], tokens[2],
                        tokens[3], tokens[4], tokens[5])
        molecule._atoms.append(atom)
    elif section == "BOND":
        bond = Mol2Bond(tokens[0], tokens[1], tokens[2], tokens[3])
        molecule._bonds.append(bond)


# endregion

# region SDF

class SdfAtom(object):
    def __init__(self, serial, x, y, z, type_):
        self.serial = serial
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.atom_type = type_

    def __str__(self):
        return "{:>6.3f} {:>6.3f} {:>6.3f} {:>2}".format(
            self.x, self.y, self.z, self.atom_type)


class SdfBond(object):
    def __init__(self, source, target, type_):
        self.origin = int(source)
        self.target = int(target)
        self.type = type_

    def __str__(self):
        return "{:>3} {:>3} {:>2}".format(
            self.origin, self.target, self.type)


class SdfMolecules(_AtomCollection, _Writable):
    """
    https://en.wikipedia.org/wiki/Chemical_table_file

    Beware that bonds use indexes starting from 1.
    """

    def __init__(self):
        _AtomCollection.__init__(self)
        self._name = None
        self._bonds = []

    def __str__(self):
        return "ATOMS\n " + \
               "\n ".join(map(str, self._atoms)) + \
               "\nBONDS\n " + \
               "\n ".join(map(str, self._bonds)) + \
               "\n"

    def write_to_stream(self, stream):
        stream.write("{}\n".format(str(self._name)))
        stream.write("  custom\n")
        stream.write("\n")
        stream.write(" {} {}  0  0  0  0              1 V2000\n".format(
            len(self._atoms), len(self._bonds)))
        for atom in self._atoms:
            stream.write(" {:>9.4f} {:>9.4f} {:>9.4f} {:<2}  0  0  0  0  0  0\n"
                         .format(atom.x, atom.y, atom.z, atom.atom_type))
        for bond in self._bonds:
            "  1  2  1  0  0  0"
            stream.write("{:>3}{:>3}  {}  0  0  0\n".format(
                bond.origin, bond.target, bond.type))
        stream.write("M  END\n")
        stream.write("$$$$\n")

    def bonds(self):
        return self._bonds

    def name(self):
        return self._name


def load_sdf_file(path):
    target = SdfMolecules()
    with open(path) as input_stream:
        # Read name
        for line in input_stream:
            name = line.rstrip()
            if len(name) > 1:
                target._name = name
                break

        # Line with source is optional, it could been an empty line.
        line = next(input_stream)
        if not line.strip() == "":
            next(input_stream)

        templine = next(input_stream)
        count_line = [templine[:3], templine[3:6]]
        atom_count = int(count_line[0])
        bond_count = int(count_line[1])

        # Read atoms.
        for index in range(atom_count):
            line = next(input_stream)
            tokens = [line[:10].strip(), line[10:20].strip(), line[20:30].strip(), line[31:33].strip()]
            target._atoms.append(SdfAtom(
                len(target._atoms) + 1,
                tokens[0], tokens[1], tokens[2], tokens[3]))

        # Read bonds.
        for index in range(bond_count):
            line = next(input_stream)
            tokens = [line[:3].strip(), line[3:6].strip(), line[6:9].strip()]
            target._bonds.append(SdfBond(
                int(tokens[0]), int(tokens[1]), tokens[2]))

        # Read rest of the file.
        for line in input_stream:
            if line == "$$$$":
                break
    return target


def load_sdf_file_property(path, property_):
    expected_line = "> <" + property_ + ">"

    with open(path) as input_stream:
        for line in input_stream:
            if line.rstrip() == expected_line:
                return next(input_stream).rstrip()

    raise RuntimeError("Missing property: " + property_ + " in: " + path)

# endregion

# region fasta


class Fasta(object):
    def __init__(self, defline):
        self.defline = defline
        self.value = ""

    def __str__(self):
        return ">{}\n".format(self.defline) + self.value


def load_fasta_file(path: str) -> typing.List[Fasta]:
    output = []
    with open(path) as input_stream:
        current_record = None
        for line in input_stream:
            line = line.rstrip()
            if line.startswith(">"):
                current_record = Fasta(line[1:].rstrip())
                output.append(current_record)
            else:
                current_record.value += line
    return output


# endregion


def add_sdf_bonds_to_pdbqt(sdf: SdfMolecules, pdbqt: PdbqtObject):
    """
    Add bonds from SDF file into PDBQT file.
    Use positions to match the atoms.
    """
    sdf_to_pdbqt_atoms = _create_atom_mapping(sdf.atoms(), pdbqt.atoms())
    # Check that we have mapping for all pdbqt atoms.
    mapping_size = len(set(sdf_to_pdbqt_atoms.values()))
    if not mapping_size == len(pdbqt.atoms()):
        raise Exception("Can't map atoms, mapped {} ouf ot {}!".format(
            mapping_size, len(pdbqt.atoms())))
    pdbqt._bonds = []
    for bond in sdf.bonds():
        origin = sdf_to_pdbqt_atoms.get(bond.origin, None)
        target = sdf_to_pdbqt_atoms.get(bond.target, None)
        if origin is None or target is None:
            # Pdbqt may not contains all molecules (hydrogens).
            continue
        pdbqt._bonds.append(PdbqtBond(origin, target))


def _create_atom_mapping(left, right, same_threshold=0.001):
    """
    Using atom indexes as stored in files. To reference an atom substract one.
    """
    mapping = {}
    for left_index, left_atom in enumerate(left):
        for right_index, right_atom in enumerate(right):
            if _same_atom_by_position(left_atom, right_atom, same_threshold):
                mapping[left_index + 1] = right_index + 1
                break
    return mapping


def _same_atom_by_position(left_atom, right_atom, same_threshold):
    max_distance = max(abs(left_atom.x - right_atom.x),
                       abs(left_atom.y - right_atom.y),
                       abs(left_atom.z - right_atom.z))
    return max_distance < same_threshold


def iter_json_lines_molecules(file_path: str):
    for model in iter_json_lines(file_path):
        for mol in model:
            yield mol


def iter_json_lines(file_path: str):
    with open(file_path, "r") as input_stream:
        for line in input_stream:
            yield json.loads(line)


def update_sdf_positions_from_pdbqt(
        sdf: SdfMolecules,
        pdbqt_mapping: PdbqtObject,
        pdbqt_target: PdbqtObject):
    """
    Use pdbqt_template to map molecules form SDF to positions in
    PDBQT mapping, then use this positions to update coordinates in
    SDF file to those in PDBQT target.
    """
    mapping = _create_atom_mapping(sdf.atoms(), pdbqt_mapping.atoms())
    sdf_atoms = sdf.atoms()
    target_atoms = pdbqt_target.atoms()
    for sdf_index, pdbqt_index in mapping.items():
        sdf_atom = sdf_atoms[sdf_index - 1]
        target_atom = target_atoms[pdbqt_index - 1]
        _transfer_position(target_atom, sdf_atom)


def _transfer_position(from_atom, to_atom):
    to_atom.x = from_atom.x
    to_atom.y = from_atom.y
    to_atom.z = from_atom.z


def remove_atoms_by_types(molecule, atom_types):
    # Get list to remove start from the biggest indexes.
    index_of_atom_to_remove = sorted([
        index + 1 for index, atom in enumerate(molecule.atoms())
        if atom.atom_type in atom_types
    ], reverse=True)
    # Delete atoms.
    atoms = molecule.atoms()
    for index_to_remove in index_of_atom_to_remove:
        atom_index = index_to_remove - 1
        del atoms[atom_index]
    # Update connections.
    bonds = molecule.bonds()
    for index in range(len(bonds) - 1, 0, -1):
        bond = bonds[index]
        if bond.origin in index_of_atom_to_remove or \
                bond.target in index_of_atom_to_remove:
            del bonds[index]
            continue
        # Update bonds indexes.
        origin_change = sum([1 for index in index_of_atom_to_remove
                             if index < bond.origin])
        target_change = sum([1 for index in index_of_atom_to_remove
                             if index < bond.target])
        bond.origin -= origin_change
        bond.target -= target_change
