#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Requires downloaded BRENDA database.
"""

import os
import logging
import collections

from libs.configuration import ENV

logger = logging.getLogger(__name__)


class BrendaEntry(object):
    def __init__(self):
        self.uniprot_ac = None
        self.substrates = []

    def __str__(self):
        return str(self.uniprot_ac) + " : " + ", ".join(self.substrates)


class _BrendaFileLoader(object):

    def __init__(self, ec_number):
        self.path = ENV["brenda_dump_path"]
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
        self.protein_synonyms = []

    def read_for_entry(self):
        with open(self.path, encoding="utf8") as in_stream:
            last_line = next(in_stream)
            for line in in_stream:
                line = line.rstrip()
                if line.startswith("\t"):
                    last_line += line[1:]
                else:
                    self._handle_line(last_line)
                    last_line = line
            self._handle_line(last_line)
        self._assign_uniprot_ac_using_synonyms()
        return list(self.output.values())

    def _handle_line(self, line):
        if line.startswith("ID"):
            self._on_id(line)
        if self.skip:
            return
        #
        line = self._remove_commentaries(line)
        line = self._remove_citations(line)
        line = self._remove_special_information(line)
        refs, line = self._extract_protein_reference(line)

        if line.startswith("PR"):
            self._on_pr(refs, line)
        elif line.startswith("SP"):
            self._on_sp(refs, line)
        elif line.startswith("SY"):
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
        reference = []
        output_line = ""
        opened = False
        buffer = ""
        for character in line:
            if character == "#":
                if opened:
                    reference.append(buffer)
                    buffer = ""
                opened = not opened
                continue
            if not opened:
                output_line += character
                continue
            if character == ",":
                reference.append(buffer)
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
        tokens = [t for t in line[3:].split(" ") if not t == ""]

        if "UniProt" not in tokens:
            name = " ".join(tokens)
            for ref in references:
                self.id_to_name[ref] = name
            return

        # TODO We can also use SwissProt
        # Example ['30'] > PR	 Sphingomonas paucimobilis P51698 SwissProt
        ref_index = tokens.index("UniProt") - 1

        uniprot_ac = tokens[ref_index]
        name = " ".join(tokens[:ref_index])
        for ref in references:
            self.output[ref].uniprot_ac = uniprot_ac
            self.id_to_name[ref] = name

    def _on_sp(self, references, line):
        if "=" not in line:
            return
        reactants = line[3:line.index("=")]
        substrates = [token.strip().rstrip()
                      for token in reactants.split(" + ")]
        for ref in references:
            self.output[ref].substrates.extend(substrates)

    def _on_sy(self, references, line):
        self.protein_synonyms.append(references)

    def _assign_uniprot_ac_using_synonyms(self):
        for key, entry in self.output.items():
            if entry.uniprot_ac is not None:
                continue
            synonym = self._find_synonym(key)
            if synonym is None:
                continue
            entry.uniprot_ac = self.output[synonym].uniprot_ac
            logger.debug("Link %s -> %s on %s ", key, synonym, self.id_to_name[key])

    def _find_synonym(self, key):
        key_name = self.id_to_name.get(key, None)

        if key_name is None:
            logger.warning("Missing name for key: %s", key)
            return None

        for group in self.protein_synonyms:

            if key not in group or len(group) == 1:
                continue

            for synonym in group:
                if synonym == key:
                    continue

                # The candidate synonym must have uniprot_ac.
                if synonym not in self.output or \
                        self.output[synonym].uniprot_ac is None:
                    continue

                synonym_name = self.id_to_name.get(synonym, None)
                if synonym_name is None:
                    logger.warning("Missing name for key: %s", synonym)
                    continue
                if key_name == synonym_name:
                    return synonym
        return None


def load_brenda_entries(ec_number):
    loader = _BrendaFileLoader(ec_number)
    return loader.read_for_entry()

# load_brenda_entries("3.8.1.5")
