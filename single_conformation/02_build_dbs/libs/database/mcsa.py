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
import pandas as pd
from lxml.cssselect import CSSSelector

from libs.utils import UrlCache
from libs.configuration import ENV


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


class McsaEntry(object):
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
    entries = [e for e in entries if (e.mcsa_id and e.enzyme_name and e.uniprot_ac and e.ec_number and e.pdb)]

    return entries


def get_mcsa_entry(mcsa_id):
    entry = McsaEntry(mcsa_id)
    mcsa_pdb = get_pdb_from_db(mcsa_id)
    chebis = _get_chebis_from_db(mcsa_pdb.lower())
    entry.reactants_chebi.extend(chebis)
    entry.catalytic_residues.extend(_get_active_residues_from_db(mcsa_pdb.lower()))
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


def query_for_ec_number_db(ec_number):
    mcsa_ids = get_mcsa_ids_from_db(ec_number)
    entries = []
    for mcsa_id in mcsa_ids:
        entry = McsaSearchEntry()
        entry.mcsa_id = mcsa_id
        entry.uniprot_ac = get_uniprot_from_db(mcsa_id)
        entry.pdb = get_pdb_from_db(mcsa_id)
        entries.append(entry)
    return entries


def get_mcsa_ids_from_db(ec_number: str) -> list:
    path = ENV["mcsa_curated_data"]
    curated_data = pd.read_csv(path)
    mcsa_entries = list(curated_data[curated_data.PDB.notnull()][curated_data.EC == ec_number]['M-CSA ID'].unique())
    return mcsa_entries


def get_uniprot_from_db(mcsa_id: str) -> str:
    path = ENV["mcsa_curated_data"]
    curated_data = pd.read_csv(path)
    uniprot_entries = curated_data[curated_data['M-CSA ID'] == mcsa_id]['Uniprot IDs'].unique()
    if uniprot_entries.size > 1:
        message = "The queried M-CSA entry %s returned more than one Uniprot entry, "\
                  "please check the values in the database manually.\nThe Uniprot "\
                  "entries are: %s" %(mcsa_id, [ue for ue in uniprot_entries])
        raise RuntimeError(message)
    elif uniprot_entries.size == 0:
        message = "The queried M-CSA entry %s returned zero Uniprot entries, "\
                  "please check the values in the database manually." % mcsa_id
        raise RuntimeError(message)
    return uniprot_entries[0]


def get_pdb_from_db(mcsa_id: str) -> str:
    path = ENV["mcsa_curated_data"]
    curated_data = pd.read_csv(path)
    pdb_entries = curated_data[curated_data['M-CSA ID'] == mcsa_id].PDB.unique()
    if pdb_entries.size > 1:
        message = "The queried M-CSA entry %s returned more than one PDB entry, "\
                  "please check the values in the database manually.\nThe PDB "\
                  "entries are: %s" %(mcsa_id, [pe for pe in pdb_entries])
        raise RuntimeError(message)
    elif pdb_entries.size == 0:
        message = "The queried M-CSA entry %s returned zero PDB entries, "\
                  "please check the values in the database manually." % mcsa_id
        raise RuntimeError(message)
    return pdb_entries[0].upper()


def _get_active_residues_from_db(pdb_id: str):
    pdb_res = get_residues_from_pdb(pdb_id)
    residues = []
    for pdb in pdb_res:
        residue = CatalyticResidue()
        residue.code = pdb[2].upper()
        residue.position = int(pdb[3])
        residue.chain = pdb[1]
        residue.types.extend(_get_residue_roles(pdb[0], pdb[3]))
        residue.description = "Not supported yet"
        residues.append(residue)
    return residues


def get_residues_from_pdb(pdb_id: str) -> set:
    path = ENV["mcsa_residues_pdb"]
    mcsa_residues = pd.read_csv(path)
    residues = set()
    pdb_residues = mcsa_residues[mcsa_residues['PDB ID'] == pdb_id]
    for i, row in pdb_residues.iterrows():
        residues.add((row['PDB ID'], row['CHAIN ID'], row['RESIDUE TYPE'],
                      row['RESIDUE NUMBER'], row['CHEMICAL FUNCTION']))
    return residues


def _get_residue_roles(pdb_id: str, res_indx: int) -> list:
    path = ENV["mcsa_residues_roles"]
    mcsa_roles = pd.read_csv(path)
    residue_roles = list()
    query_res = mcsa_roles[mcsa_roles['PDB ID'] == pdb_id][mcsa_roles['RESIDUE NUMBER'] == res_indx]
    for i, row in query_res.iterrows():
        residue_roles.append(row['ROLE'])
    return residue_roles


def _get_chebis_from_db(pdb_id: str) -> list:
    path = ENV["mcsa_curated_data"]
    curated_data = pd.read_csv(path)
    chebis = list(curated_data[curated_data.PDB == pdb_id][curated_data['residue/reactant/product/cofactor'] == 'reactant']['resid/chebi id'].unique())
    return chebis


def add_substrate_from_db(uniprot_entry):
    path = ENV["mcsa_db"]
    curated_data = pd.read_csv(path)
    uniprot_e = curated_data[curated_data['Uniprot IDs'] == uniprot_entry.ac()]
    reactants = uniprot_e[uniprot_e['residue/reactant/product/cofactor'] == 'reactant']
    substrates = reactants['resid/chebi id'].unique()
    for s in substrates:
        if s != '15377':
            uniprot_entry.prop("substrate").add(s).source("m-csa")