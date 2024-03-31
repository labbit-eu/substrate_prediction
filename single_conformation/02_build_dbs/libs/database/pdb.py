#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Protein Data Bank database
"""

import urllib.request
import xml.sax

from libs import utils


class PdbDescription(object):

    def __init__(self, structure_id):
        self.structure_id = structure_id
        self.title = None
        self.resolution = None

    def __str__(self):
        return self.structure_id + " (" + self.resolution + ") : " + self.title


class PdbDescriptionFactory(xml.sax.ContentHandler):
    pdb_rest_url = "https://www.rcsb.org/pdb/rest/"
    cache_name = "pdb-description"

    def __init__(self):
        super().__init__()
        self.value = None

    def create(self, structure_id):
        self.value = PdbDescription(structure_id)

        url = self.pdb_rest_url + "describePDB?structureId=%s" % structure_id
        with urllib.request.urlopen(url) as stream:
            xml.sax.parse(stream, self)

        return self.value

    @staticmethod
    def to_json(data):
        return data.__dict__

    @staticmethod
    def from_json(data):
        value = PdbDescription(data["structure_id"])
        value.title = data["title"]
        value.resolution = data["resolution"]
        return value

    def startElement(self, tag, attributes):
        if not tag == "PDB":
            return
        self.value.title = attributes["title"]
        self.value.resolution = attributes["resolution"]


class PdbLigand(object):

    def __init__(self):
        self.chemical_id = None
        self.chemical_name = None
        self.type = None
        self.formula = None
        self.smiles = None

    def __str__(self):
        return self.chemical_name + " : " + self.smiles


class PdbLigandsFactory(xml.sax.ContentHandler):
    pdb_rest_url = "https://www.rcsb.org/pdb/rest/"
    cache_name = "pdb-ligand"

    def __init__(self):
        self.open_value = None
        self.values = []
        self.xml_open_tags = []

    def create(self, structure_id):
        self.values = []

        url = self.pdb_rest_url + "ligandInfo?structureId=%s" % structure_id
        with urllib.request.urlopen(url) as stream:
            xml.sax.parse(stream, self)

        return self.values

    @staticmethod
    def to_json(data):
        return [item.__dict__ for item in data]

    @staticmethod
    def from_json(data):
        values = []
        for item in data:
            value = PdbLigand()
            value.chemical_id = item["chemical_id"]
            value.chemical_name = item["chemical_name"]
            value.type = item["type"]
            value.formula = item["formula"]
            value.smiles = item["smiles"]
            values.append(value)
        return values

    def startElement(self, tag, attributes):
        self.xml_open_tags.append(tag)
        if not tag == "ligand":
            return
        self.open_value = PdbLigand()
        self.open_value.chemical_id = attributes["chemicalID"]
        self.open_value.type = attributes["type"]

    def endElement(self, tag):
        self.xml_open_tags.pop()
        if tag == "ligand":
            self.values.append(self.open_value)
            self.open_value = None

    def characters(self, content):
        if self.open_value is None:
            return
        open_tag = self.xml_open_tags[-1]
        if open_tag == "formula":
            self.open_value.formula = content
        elif open_tag == "smiles":
            self.open_value.smiles = content
        elif open_tag == "chemicalName":
            self.open_value.chemical_name = content


def download_pdb(pdb_id):
    url = "https://files.rcsb.org/download/" + pdb_id + ".pdb"
    return utils.fetch_url(pdb_id, url, "pdb-files")


def fetch_fasta(pdb_id):
    url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do" + \
          "?structureIdList=" + pdb_id + \
          "&compressionType=uncompressed"
    return utils.fetch_url(pdb_id, url, "pdb-fasta")


def is_fasta_ok(path):
    with open(path, 'r') as fasta_stream:
        first_line = fasta_stream.readline()
        return first_line.startswith('>')

