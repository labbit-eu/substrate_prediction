#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Store definition of all data objects.
"""

import os
import json
import typing
import logging

logger = logging.getLogger(__name__)


# region Base objects definitions

class ProvenanceAware(object):
    """
    Add options to store provenance data.
    """

    def __init__(self):
        self._provenance = []

    def source(self, source: str):
        self.provenance(source)

    def provenance(self, source: str, additional: object = None):
        self._provenance.append({
            "source": source,
            "additional": additional
        })

    def sources(self):
        return [item["source"] for item in self._provenance]


class Value(ProvenanceAware):
    """
    Represent a value.
    """

    def __init__(self, value):
        super(Value, self).__init__()
        self._value = value

    def value(self):
        return self._value

    def replace(self, value) -> "Value":
        self._value = value
        return self


class Property(object):
    """
    Collection of values under given name.
    """

    def __init__(self):
        self._values = []

    def set(self, value) -> Value:
        new_value = Value(value)
        self._values = [new_value]
        return new_value

    def add(self, value, ignore_none=True) -> Value:
        if value is None and ignore_none:
            # Return value, which is of no use but user does not have to
            # check.
            return Value(None)
        new_value = Value(value)
        self._values.append(new_value)
        return new_value

    def add_array(self, values, ignore_none=True):
        for value in values:
            self.add(value, ignore_none)

    def get(self) -> typing.List[Value]:
        return self._values

    def first(self) -> typing.Optional[Value]:
        if len(self._values) == 0:
            return None
        else:
            return self._values[0]

    def no_value(self):
        return len(self._values) == 0

    def filter(self, predicate):
        self._values = [value for value in self._values if predicate(value)]

    def as_array(self):
        return [val.value() for val in self._values]

    def map(self, transformer):
        for value in self._values:
            value._value = transformer(value._value)


class Entry(ProvenanceAware):
    """
    An object with values stores in properties.
    """

    def __init__(self, user_type, user_id):
        super(Entry, self).__init__()
        self._ref_id = None  # Is set as added into a database.
        self._user_type = user_type
        self._user_id = user_id
        self._properties = {}

    def __str__(self):
        output = "Entry: {} {} [{}]\n".format(
            self._user_type, self._user_id, ", ".join(self.sources()))

        for name, prop in self._properties.items():
            output += "  {}\n".format(name)
            for value in prop.get():
                output += "    {} [{}]\n".format(
                    str(value.value()), ", ".join(value.sources()))

        return output

    def _init_params(self):
        """
        :return: Parameters passed to __init__ method of this instance.
        """
        return [self._user_type, self._user_id]

    def ref(self):
        """
        :return: Value that can be used as a reference to this object.
        """
        return self._ref_id

    @staticmethod
    def type() -> str:
        raise NotImplementedError()

    def prop(self, name: str) -> Property:
        if name not in self._properties:
            self._properties[name] = Property()
        return self._properties[name]

    def peek_prop(self, name: str) -> typing.Optional[Property]:
        return self._properties.get(name, None)

    def properties(self):
        return self._properties


# endregion


class UniprotEntry(Entry):

    def __init__(self, ac: str):
        super(UniprotEntry, self).__init__(self.type(), ac)
        self._ac = ac

    def _init_params(self):
        return [self._ac]

    @staticmethod
    def type():
        return "uniprot"

    def ac(self):
        return self._ac

    def ec_number(self):
        return self.prop("ec_number")

    def name(self):
        return self.prop("name")

    def pdb(self):
        return self.prop("pdb")

    def chembl(self):
        return self.prop("chembl")

    def residues(self):
        return self.prop("residues")

    def substrate(self):
        return self.prop("substrate")

    def mcsa(self):
        return self.prop("m-csa")


class McsaEntry(Entry):

    def __init__(self, mcsa_id: str, chain=None):
        super(McsaEntry, self).__init__(self.type(), mcsa_id)
        self._mcsa_id = mcsa_id
        self._active_chain = chain

    def _init_params(self):
        return [self._mcsa_id]

    @staticmethod
    def type():
        return "m-csa"

    def mcsa_id(self):
        return self._mcsa_id

    def active_chain(self):
        return self._active_chain

    def pdb(self):
        return self.prop("pdb")

    def residues(self):
        return self.prop("residues")

    def chebi(self):
        """Fragment definition."""
        return self.prop("chebi")


class SubstrateEntry(Entry):

    def __init__(self, database: str, id_in_database: str):
        user_id = database + "-" + id_in_database
        super(SubstrateEntry, self).__init__(self.type(), user_id)
        # self.database().add(database)

    def _init_params(self):
        return self._user_id.split("-")

    @staticmethod
    def type():
        return "substrate"

    def database(self):
        return self.prop("database")

    def name(self):
        return self.prop("name")

    def query(self):
        return self.prop("query")

    def has_no_pdbqt(self):
        """
        Set and true if we fail to create a PDBQT structure.
        """
        return self.prop("has_no_pdbqt")


class MappingEntry(Entry):

    def __init__(self, source_pdb: str, source_chain: str,
                 target_pdb: str, target_chain: str, mcsa: str):
        user_id = "-".join([source_pdb, source_chain,
                            target_pdb, target_chain,
                            mcsa])
        super(MappingEntry, self).__init__(self.type(), user_id)

    def _init_params(self):
        return self._user_id.split("-")

    @staticmethod
    def type():
        return "mapping"

    def mcsa(self):
        return self.prop("m-csa")

    def source_pdb(self):
        return self.prop("source_pdb")

    def source_chain(self):
        return self.prop("source_chain")

    def target_pdb(self):
        return self.prop("target_pdb")

    def target_chain(self):
        return self.prop("target_chain")

    def is_valid(self):
        """
        False is the mapping is invalid.
        """
        return self.prop("is_valid")

    def has_no_pdbqt(self):
        """
        Set and true if we fail to create a PDBQT structure.
        """
        return self.prop("has_no_pdbqt")

    def offset(self):
        return self.prop("offset")

    def fasta_to_mcsa(self):
        return self.prop("fasta_to_mcsa")


class ComplexEntry(Entry):

    def __init__(self, pdb: str, chain: str, substrate_ref: str):
        self._ctor_params = [pdb, chain, substrate_ref]
        user_id = "-".join([pdb, chain, substrate_ref])
        super(ComplexEntry, self).__init__(self.type(), user_id)

    def _init_params(self):
        return self._ctor_params

    @staticmethod
    def type():
        return "complex"

    def substrate(self):
        return self.prop("substrate")

    def protein(self):
        return self.prop("protein")

    def mapping(self):
        return self.prop("mapping")


class Fragment(object):
    class Atom(object):
        def __init__(self):
            self.x = 0
            self.y = 0
            self.z = 0
            self.serial = None
            self.name = None

    def __init__(self):
        self.mapping = None
        self.substrate = None
        self.name = None
        self.subdir_name = None
        self.model_index = None
        self.fragment_index = None
        self.vina_score = None
        self.atoms = None
        # Computed
        self.coulomb_score = None
        self.distances = None
        self.accessibility = None
        # Other temporary.
        self.props = {}
        

# region Database

class Database(object):

    def __init__(self):
        self._entries = {}

    def add(self, entry: Entry):
        ref_key = self.create_ref_id(entry.type(), entry._user_id)
        entry._ref_id = ref_key
        self._entries[ref_key] = entry

    @staticmethod
    def create_ref_id(type_: str, user_id: str):
        return "@" + type_ + "," + user_id

    def get(self, type_: str, user_id: str) -> Entry:
        key = self.create_ref_id(type_, user_id)
        return self._entries.get(key, None)

    def get_or_add(self, entry: Entry) -> Entry:
        ref_key = self.create_ref_id(entry.type(), entry._user_id)
        if ref_key in self._entries:
            return self._entries[ref_key]
        #
        self._entries[ref_key] = entry
        entry._ref_id = ref_key
        return entry

    def iter(self, type_to_iter=None) -> typing.Generator[Entry, None, None]:
        try:
            type_to_iter = type_to_iter.type()
        except AttributeError:
            pass

        for entry in self._entries.values():
            if type_to_iter is None or entry.type() == type_to_iter:
                yield entry

    def ref(self, id_: str) -> Entry:
        return self._entries.get(id_, None)


class ReadOnlyDatabase(object):

    def __init__(self, directory):
        self._directory = directory

    def read(self, entry_type, ref_id: str) -> Entry:
        path = os.path.join(
            self._directory, entry_type.type(), ref_id + ".json")
        with open(path) as in_stream:
            obj_json = json.load(in_stream)
        return load_entry_from_json(entry_type, obj_json)

    def iter(self, entry_type) -> typing.Generator[Entry, None, None]:
        type_dir = os.path.join(self._directory, entry_type.type())
        for file_name in os.listdir(type_dir):
            path = os.path.join(type_dir, file_name)
            with open(path) as in_stream:
                obj_json = json.load(in_stream)
            yield load_entry_from_json(entry_type, obj_json)


def create_in_memory() -> Database:
    return Database()


def create_read_only_from_directory(directory) -> ReadOnlyDatabase:
    return ReadOnlyDatabase(directory)


def save_database(directory: str, database: Database):
    logger.debug("Saving database to: %s", directory)
    os.makedirs(directory, exist_ok=True)
    for item in database.iter():
        user_type = item._user_type
        ref_id = item._ref_id
        os.makedirs(os.path.join(directory, user_type), exist_ok=True)
        path = os.path.join(directory, user_type, ref_id + ".json")
        with open(path, "w") as out_stream:
            dump_entry(item, out_stream, indent=2)


def is_ref(value):
    return value.value().startswith("@")


def filter_prop(entries, property_name, property_value):
    output = []
    for entry in entries:
        prop = entry.peek_prop(property_name)
        if prop is None:
            continue
        for value in prop.get():
            if value.value() == property_value:
                output.append(entry)
                break
    return output


# endregion

# region Input/Output

def dump_entry(entry, stream, indent=None):
    return json.dump(entry, stream, indent=indent, cls=_EntryJsonEncoder)


class _EntryJsonEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Entry):
            return self._on_entry(o)
        elif isinstance(o, Property):
            return self._on_property(o)
        elif isinstance(o, Value):
            return self._on_value(o)
        elif isinstance(o, Fragment):
            return self._on_fragment(o)
        elif isinstance(o, Fragment.Atom):
            return self._on_fragment_atom(o)
        else:
            return json.JSONEncoder.default(self, o)

    @staticmethod
    def _on_entry(o: Entry):
        return {
            "__init__": o._init_params(),
            "ref_id": o._ref_id,
            "properties": o._properties,
            "provenance": _EntryJsonEncoder._provenance(o._provenance)
        }

    @staticmethod
    def _on_property(o: Property):
        return o._values

    @staticmethod
    def _on_value(o: Value):
        return {
            "value": o._value,
            "provenance": _EntryJsonEncoder._provenance(o._provenance)
        }

    @staticmethod
    def _provenance(provenance):
        output = []
        for item in provenance:
            if item["additional"] is None:
                output.append({
                    "source": item["source"]
                })
            else:
                output.append(item)
        return output

    @staticmethod
    def _on_fragment(o: Fragment):
        return {
            "name": o.name,
            "atoms": o.atoms,
            "ref": {
                "substrate": o.substrate,
                "mapping": o.mapping
            },
            "model_index": o.model_index,
            "fragment_index": o.fragment_index,
            "subdir_name": o.subdir_name,
            "vina_score": o.vina_score,
            # Computed
            "coulomb_score": o.coulomb_score,
            "distances": o.distances,
            # Temporary
            "props": o.props
        }

    @staticmethod
    def _on_fragment_atom(o: Fragment.Atom):
        return {
            "x": o.x,
            "y": o.y,
            "z": o.z,
            "serial": o.serial,
            "name": o.name
        }


def load_entry_from_json(entry_type, json_object):
    output = entry_type(*json_object["__init__"])
    output._ref_id = json_object["ref_id"]
    _load_provenance(json_object["provenance"], output)
    for key, values in json_object["properties"].items():
        prop = output.prop(key)
        for value in values:
            instance = prop.add(value["value"])
            _load_provenance(value["provenance"], instance)
    return output


def _load_provenance(provenance, obj: ProvenanceAware):
    for item in provenance:
        obj.provenance(item["source"], item.get("additional", None))

# endregion
