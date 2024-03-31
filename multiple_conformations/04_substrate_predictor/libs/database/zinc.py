#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Inspect:
* http://files.docking.org/

There are multiple download options see (mol2, solvation)
http://zinc15.docking.org/substances/ZINC000000001748/


"""

import urllib

import lxml
import lxml.html
from lxml.cssselect import CSSSelector

from utils import UrlCache


class ZincMolecule(object):
    def __init__(self, zinc_id):
        self.zinc_id = zinc_id
        self.mol2_links = []
        self.chebi_id = None
        self.chembl_id = None
        self.smiles = None

    def __str__(self):
        return self.zinc_id + " " + \
               str(self.chebi_id).rjust(15) + " " + \
               str(self.chembl_id).rjust(15) + " " + \
               self.smiles


def query_for_zinc_id(query):
    query = urllib.parse.urlencode({
        "q": query
    })
    url = "http://zinc15.docking.org/substances/search/?" + query
    item_selector = CSSSelector("h4.zinc-id > a")

    # TODO This is not save, we should replace UrlCache.
    name = _encode_smiles_for_hdd(query)
    with UrlCache(name, url, "zinc-search") as request:
        document = lxml.html.parse(request)

    return [_link_to_zinc_id(item.get("href"))
            for item in item_selector(document)]


def _encode_smiles_for_hdd(smiles):
    return smiles


def _link_to_zinc_id(link):
    return link[12:-1]


# TODO Convert to factory.
def fetch_zinc_molecule(zinc_id):
    url = "http://zinc15.docking.org/substances/" + zinc_id + "/"
    with UrlCache(zinc_id, url, "zinc-detail") as request:
        document = lxml.html.parse(request)

    molecule = ZincMolecule(zinc_id)

    download_selector = CSSSelector("td[title=Downloads] ul a")
    for element in download_selector(document):
        link = element.get("href")
        if link.endswith(".mol2.gz"):
            molecule.mol2_links.append(link)

    smiles_selector = CSSSelector("#substance-smiles-field")
    molecule.smiles = smiles_selector(document)[0].get("value")

    properties_selector = CSSSelector(
        ".annotated .catalog-items-map ul.list-inline a")
    for element in properties_selector(document):
        link = element.get("href")
        if "CHEBI" in link:
            molecule.chebi_id = _link_to_chebi_id(link)
        elif "CHEMBL" in link:
            # TODO There are links for multiple versions.
            molecule.chembl_id = _link_to_chembl_id(link)

    return molecule


def _link_to_chebi_id(link):
    return link[link.rfind("=") + 1:]


def _link_to_chembl_id(link):
    return link[link.rfind("/") + 1:]
