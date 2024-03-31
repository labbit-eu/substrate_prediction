#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mechanism and Catalytic Site Atlas:
    http://www.ebi.ac.uk/thornton-srv/m-csa/
Data source URL:
    http://www.ebi.ac.uk/thornton-srv/m-csa/search/?s={ec number}
"""

import lxml
import lxml.html
from lxml.cssselect import CSSSelector

from libs.utils import UrlCache


class McsaSearchEntry(object):

    def __init__(self):
        self.mcsa_id = None
        self.enzyme_name = None
        self.uniprot_ac = None
        self.ec_number = None
        self.pdb = None

    def __str__(self):
        return str(self.mcsa_id).rjust(5) + \
               str(self.uniprot_ac).rjust(15) + \
               str(self.pdb).rjust(10) + \
               str(self.enzyme_name).rjust(44)


class CatalyticResidue(object):
    def __init__(self):
        self.code = None
        self.position = None
        self.types = []
        self.chain = None

    def __str__(self):
        return self.code + " " + \
               str(self.position).rjust(3) + " " + \
               str(self.chain) + " : " + \
               ", ".join(self.types)


class ContainerMcsaEntry(object):
    """
    Includes only additional information to McsaSearchEntry.
    """

    def __init__(self, mcsa_id):
        self.mcsa_id = mcsa_id
        self.reactants_chebi = []
        self.catalytic_residues = []

    def __str__(self):
        return self.mcsa_id + " : " + " ".join(self.reactants_chebi) + \
               " ( " + \
               " | ".join([str(res) for res in self.catalytic_residues]) + \
               " )"


def query_for_ec_number(ec_number):
    url = "http://www.ebi.ac.uk/thornton-srv/m-csa/search/?s=" + ec_number
    with UrlCache(ec_number, url, "mcsa-search") as stream:
        document = lxml.html.parse(stream)

    row_selector = CSSSelector('.table-container tbody tr')
    cell_selector = CSSSelector('td')
    entries = []
    for row in row_selector(document):
        entry = McsaSearchEntry()
        for tag in cell_selector(row):
            tag_class = tag.get("class")
            tag_value = tag.text_content()
            if tag_class == "macie_id":
                entry.mcsa_id = tag_value
            elif tag_class == "enzyme_name":
                entry.enzyme_name = tag_value
            elif tag_class == "protein":
                entry.uniprot_ac = tag_value
            elif tag_class == "ec":
                entry.ec_number = tag_value
            elif tag_class == "pdb":
                entry.pdb = tag_value.upper()
        entries.append(entry)

    return entries


def get_mcsa_entry(mcsa_id):
    url = "http://www.ebi.ac.uk/thornton-srv/m-csa/entry/" + mcsa_id + "/"
    with UrlCache(mcsa_id, url, "mcsa-entry") as stream:
        document = lxml.html.parse(stream)

    entry = ContainerMcsaEntry(mcsa_id)
    reactants_selector = CSSSelector("#reaction-boxes div")
    chebi_link_selector = CSSSelector("div.compound-chebi a")
    for element in reactants_selector(document):
        if "compound-box" in element.get("class"):
            for link_element in chebi_link_selector(element):
                href = link_element.get("href")
                if "chebiId" in href:
                    entry.reactants_chebi.append(_link_to_chebi_id(href))
                    break
        elif "symbol-box" in element.get("class"):
            symbol_class = element[0].get("class")
            if "arrow" in symbol_class:
                # From now on there are products.
                break

    residues_selector = CSSSelector("#mech1summary tr")
    for element in residues_selector(document):
        entry.catalytic_residues.extend(_parse_residues_information(element))

    return entry


def _link_to_chebi_id(link):
    return link[link.rfind(":") + 1:]


def _parse_residues_information(element):
    names = element[0].text_content()
    description = element[1].text_content()
    types = element[2].text_content().strip().split(", ")

    output = []
    for name in names.split(", "):
        residue = CatalyticResidue()
        residue.code = name[0:3].upper()
        residue.position = int(name[3:-2])
        residue.chain = name.rstrip()[-1].upper()
        residue.description = description
        residue.types = set(types)
        output.append(residue)
    return output
