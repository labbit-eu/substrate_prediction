#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Requires downloaded BRENDA database.
"""

import os
import typing
import logging
import collections
from libs.database import uniprot as uniprot_db

logger = logging.getLogger(__name__)


class BrendaEntry(object):
    def __init__(self):
        self.brenda_id = None
        self.species = None
        self.recommended_name = None
        self.systematic_name = None
        self.uniprot_ac = []
        self.prot_names = []
        self.substrates = []
        self.substrates_smiles = {}

    def __str__(self):
        return str(self.brenda_id) + " : " + self.species + " : " +\
               self.recommended_name + " : " + self.systematic_name + " : " +\
               ", ".join(self.prot_names) + " : " + ", ".join(self.uniprot_ac)+\
               " : " + ", ".join(self.substrates)


class _BrendaFileLoader(object):

    def __init__(self, ec_number):
        self.path = "/HDD4TB/Carlos/Projects/Project_binding_modes/filter_binding_modes/db/brenda_download.txt"
        if self.path is None:
            raise RuntimeError("Invalid configuration missing "
                               "`env.brenda_dump_path` property.")
        if not os.path.exists(self.path):
            message = "Missing BRENDA database file.\n" \
                      "Please download it from " \
                      "http://www.brenda-enzymes.org/" \
                      "download_brenda_without_registration.php " \
                      "and save it as `{}`".format(self.path)
            raise RuntimeError(message)
        self.ec_number = ec_number
        self.skip = True
        self.output = collections.defaultdict(BrendaEntry)
        self.id_to_name = {}
        self.recommended_name = None
        self.systematic_name = None

    def read_for_entry(self):
        with open(self.path, encoding="utf8") as in_stream:
            last_line = next(in_stream)
            for line in in_stream:
                line = line.rstrip()
                if line.startswith("\t"):
                    last_line += " " + line[1:]
                else:
                    self._handle_line(last_line)
                    last_line = line
            self._handle_line(last_line)
        return list(self.output.values())

    def _handle_line(self, line):
        if line.startswith("ID"):
            self._on_id(line)
        if self.skip:
            return
        line = self._remove_citations(line)
        line = self._remove_special_information(line)
        refs, line = self._extract_protein_reference(line)
    
        if line.startswith("PR\t"):
            self._on_pr(refs, line)
        elif line.startswith("RN\t"):
            self._on_rn(line)
        elif line.startswith("SN\t"):
            self._on_sn(line)
        elif line.startswith("SP\t"):
            self._on_sp(refs, line)
        elif line.startswith("SY\t"):
            self._on_sy(refs, line)

    @staticmethod
    def _remove_commentaries(line):
        output = ""
        skip = 0
        for character in line:
            skip += character == "("
            if skip == 0:
                output += character
            skip -= character == ")"
        return output

    @staticmethod
    def _remove_citations(line):
        output = ""
        skip = 0
        for character in line:
            skip += character == "<"
            if skip == 0:
                output += character
            skip -= character == ">"
        return output

    @staticmethod
    def _remove_special_information(line):
        output = ""
        skip = 0
        for character in line:
            skip += character == "{"
            if skip == 0:
                output += character
            skip -= character == "}"
        return output

    @staticmethod
    def _extract_protein_reference(line):
        reference = set()
        output_line = ""
        opened = False
        buffer = ""
        for character in line:
            if character == "#":
                if opened:
                    reference.add(buffer)
                    buffer = ""
                opened = not opened
                continue
            if not opened:
                output_line += character
                continue
            if character == "," or character == " ":
                reference.add(buffer)
                buffer = ""
            else:
                buffer += character

        return reference, output_line

    def _on_id(self, line):
        tokens = line[3:].split(" ")
        if not tokens[0] == self.ec_number:
            self.skip = True
            return
        if len(tokens) > 1:
            raise RuntimeError(
                "Please implement handling for BRENDA ID line" + line)
        self.skip = False

    def _on_pr(self, references, line):
        line = self._remove_commentaries(line)
        tokens = [t for t in line[3:].split(" ") if not t == ""]

        if self._get_id_index(tokens) == -1:
            name = " ".join(tokens)
            for ref in references:
                self.output[ref].brenda_id = int(ref)
                self.output[ref].species = name
            return

        ref_index = self._get_id_index(tokens)
        uniprot_ac = tokens[ref_index]
        name = " ".join(tokens[:ref_index])
        for ref in references:
            self.output[ref].brenda_id = int(ref)
            self.output[ref].uniprot_ac.append(uniprot_ac)
            self.output[ref].species = name

    def _on_rn(self, line):
        tokens = [t for t in line[3:].split(" ") if not t == ""]
        self.recommended_name = " ".join(tokens)
        for key in self.output.keys():
            self.output[key].recommended_name = self.recommended_name

    def _on_sn(self, line):
        name = line[3:].split(":")[-1]
        self.recommended_name = name
        for key in self.output.keys():
            self.output[key].systematic_name = self.recommended_name

    def _on_sp(self, references, line):
        if "=" not in line:
            return
        reactants = line[3:line.index("=")]
        substrates = [token.strip().rstrip()
                      for token in reactants.split(" + ")]
        for ref in references:
            self.output[ref].substrates.extend(substrates)

    def _on_sy(self, references, line):
        line = self._remove_commentaries(line)
        tokens = [t for t in line[3:].split(" ") if not t == ""]
        name = " ".join(tokens)
        for ref in references:
            self.output[ref].prot_names.append(name)

    def _get_id_index(self, tokens):
        for i, token in enumerate(tokens):
            if token.lower() == "uniprot" or token.lower() == "uniprot," or \
               token.lower() == "swissprot" or token.lower() == "swissprot,":
                return i-1
        return -1


def load_brenda_entries(ec_number):
    loader = _BrendaFileLoader(ec_number)
    return loader.read_for_entry()


def brenda_finder(brentries):
    for entry in brentries:
        if len(entry.substrates) > 0:
            uniprot_ac = uniprot_db.search_by_sp_and_name(entry.species,
                                                          entry.recommended_name)
            if uniprot_ac is not None:
                entry.uniprot_ac.append(uniprot_ac)
            else:
                uniprot_ac = uniprot_db.search_by_sp_and_name(entry.species,
                                                          entry.systematic_name)
                if uniprot_ac is not None:
                    entry.uniprot_ac.append(uniprot_ac)
                else:
                    for name in entry.prot_names:
                        uniprot_ac = uniprot_db.search_by_sp_and_name(entry.species,
                                                                      name)
                        if uniprot_ac is not None:
                            entry.uniprot_ac.append(uniprot_ac)
