#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of Uniprot record with links to ChEMBL and PDB
    http://www.uniprot.org/uniprot/P27652
"""

import os
import json
import http
import typing
import urllib.request
import xml.sax
import logging
from bs4 import BeautifulSoup

from libs import utils
from libs import paths
from libs.configuration import TASK

PDB_REF = "PDB"
CHEMBL_REF = "ChEMBL"
TCDB_REF = "TCDB"

logger = logging.getLogger(__name__)


class ActiveSite(object):
    """
    Representation of active site (residue).
    """

    def __init__(self, description):
        self.description = description
        self.position = None

    def __str__(self):
        return self.position + ":" + str(self.description)


class UniprotEntry(object):
    REFERENCES = [PDB_REF, CHEMBL_REF, TCDB_REF]

    def __init__(self, ac):
        self.ac = ac
        # TODO Replace "dictionary" entry with class objects.
        self.db_references = []
        self.ec_number = None
        self.active_sites = []
        self.recommended_name = None
        # TODO Extract Catalytic activity
        # http://www.uniprot.org/uniprot/B8YPW8

    def __str__(self):

        return self.ac.rjust(12) + \
               str(list(self.references_id(CHEMBL_REF))).rjust(20) + \
               str(self.ec_number).rjust(10) + \
               str([str(x) for x in self.active_sites]).rjust(65)

    def references_id(self, reference_type=None):
        if type is None:
            for entry in self.db_references:
                yield entry["id"]
        else:
            for entry in self.db_references:
                if not entry["type"] == reference_type:
                    continue
                yield entry["id"]


class UniprotEntryFactory(xml.sax.ContentHandler):
    cache_name = "uniprot"

    def __init__(self):
        self.value = None
        self.xml_open_tags = []
        self.xml_open_active = None

    def create(self, ac):
        self.value = UniprotEntry(ac)

        url = "http://www.uniprot.org/uniprot/%s.xml" % ac
        with urllib.request.urlopen(url) as stream:
            xml.sax.parse(stream, self)

        value = self.value
        self.value = None
        return value

    @staticmethod
    def to_json(data):
        return {
            "ac": data.ac,
            "recommended_name": data.recommended_name,
            "db_references": data.db_references,
            "ec_number": data.ec_number,
            "active_sites": [item.__dict__ for item in data.active_sites]
        }

    @staticmethod
    def from_json(data):
        value = UniprotEntry(data["ac"])
        value.recommended_name = data["recommended_name"]
        value.db_references = data["db_references"]
        value.ec_number = data["ec_number"]
        value.active_sites = []
        for item in data["active_sites"]:
            active_site = ActiveSite(item["description"])
            active_site.position = item["position"]
            value.active_sites.append(active_site)
        return value

    def startElement(self, tag, attributes):
        self.xml_open_tags.append(tag)
        if tag == "dbReference" and self.xml_open_tags[-2] == "entry":
            if attributes["type"] in UniprotEntry.REFERENCES:
                self.value.db_references.append({
                    "id": attributes["id"],
                    "type": attributes["type"]
                })
        elif tag == "feature" and attributes["type"] == "active site":
            self.xml_open_active = ActiveSite(attributes.get("description", None))
        elif tag == "position" and self.xml_open_active is not None:
            self.xml_open_active.position = attributes["position"]

    def endElement(self, tag):
        self.xml_open_tags.pop()
        if tag == "feature" and self.xml_open_active is not None:
            self.value.active_sites.append(self.xml_open_active)
            self.xml_open_active = None

    def characters(self, content):
        open_tag = self.xml_open_tags[-1]
        if open_tag == "ecNumber":
            if self.value.ec_number is None:
                self.value.ec_number = content
            else:
                logger.debug("More than one EC found for protein %s. "
                             "Using only first one: %s, other EC: %s" %(self.value.ac, self.value.ec_number, content))
        if open_tag == "fullName" and \
                self.xml_open_tags[-2] == "recommendedName":
            self.value.recommended_name = content


class SearchResultItem(object):
    """
    Uniprot search result item.
    """

    def __init__(self, ac):
        self.ac = ac
        self.identity = None

    def __str__(self):
        return self.ac.ljust(12) + "(" + self.identity + ")"


class UniprotSearchRemoteResultHandler(xml.sax.ContentHandler):
    """
    Parser of Uniprot search result file with automatic command:
    wget "https://web.expasy.org/cgi-bin/blast/blast.pl?seq=SEQ&format=xml&prot_db1=UniProtKB&matrix=BLOSUM62&showal=
    250&ethr=10&Filter=F&Gap=T" -O "test.xml" --timeout=6000
    """
    
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.result_item = None
        self.identity_element_open = False
        self.results = []
    
    def startElement(self, name, attributes):
        if name == "Hit":
            self.result_item = SearchResultItem("temporal")
            self.result_item.identity = 0
        elif name == "Hit_accession" and self.result_item is not None:
            self.result_item.ac = "temporal"
        elif name == "Hsp_identity" and self.result_item is not None:
            self.result_item.identity = None
    
    def characters(self, content):
        if self.result_item is not None:
            if self.result_item.ac == "temporal":
                self.result_item.ac = content.strip()
            if self.result_item.identity is None:
                self.result_item.identity = content.strip()
    
    def endElement(self, name):
        if name == "Hit":
            self.result_item.identity = str(float(self.result_item.identity)/298.0*100)
            self.results.append(self.result_item)
            self.result_item = None
    
    def __iter__(self):
        return iter(self.results)


class UniprotSearchResultHandler(xml.sax.ContentHandler):
    """
    Parser of Uniprot search result file.
    """

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.result_item = None
        self.identity_element_open = False
        self.results = []

    def startElement(self, tag, attributes):
        if tag == "hit":
            self.result_item = SearchResultItem(attributes["ac"])
        elif tag == "identity":
            self.identity_element_open = True

    def endElement(self, tag):
        if tag == "hit":
            self.results.append(self.result_item)
            self.result_item = None
        elif tag == "identity":
            self.identity_element_open = False

    def characters(self, content):
        if self.identity_element_open and self.result_item is not None:
            self.result_item.identity = content.strip()

    def __iter__(self):
        return iter(self.results)


def fetch_search_results(sequence) -> typing.List[SearchResultItem]:
    # TODO #8
    path = os.path.join(paths.manual(), "uniprot_fasta_query.xml")

    if TASK["knowledge_base"]["uniprot_search_remote"]:
        search_result = UniprotSearchRemoteResultHandler()
        if not os.path.exists(path):
            logger.info("Searching similar sequences in Uniprot database at web.expasy.org. This might take a while...")
            base_url = "https://web.expasy.org/cgi-bin/blast/blast.pl?"
            args = "seq={}&format=xml&prot_db1=UniProtKB&matrix=BLOSUM62&showal=250&ethr={}&Filter=F&Gap=T".format(
                sequence, TASK["knowledge_base"]["uniprot_search_evalue"])

            uniprot_url = base_url + args

            if not utils.blast_uniprot(uniprot_url, path):
                message = "Missing Uniprot fasta search file.\n" \
                          "Please check internet connection or set the variable " \
                          "'uniprot_search_remote' to False and proceed manually."
                raise RuntimeError(message)
    else:
        if not os.path.exists(path):
            message = "Missing Uniprot fasta search file.\n" \
                      "Please query Uniprot with the input protein sequence " \
                      "and save results as xml " \
                      "into " + path + \
                      "\n\nInput protein sequence:\n\n" + sequence
            raise RuntimeError(message)
        search_result = UniprotSearchResultHandler()
    xml.sax.parse(path, search_result)
    return search_result


def fetch_uniprot_entry(uniprot_ac) -> UniprotEntry:
    uniprot_entry_factory = UniprotEntryFactory()
    try:
        with utils.ObjectCache(uniprot_entry_factory, uniprot_ac) as uniprot_entry:
            pass
    except:
        logger.exception("Ignoring invalid (error reading file) "
                         "Uniprot record: " + uniprot_ac)
        return None

    return uniprot_entry


def search_by_sp_and_name(species_name, prot_name):
    species = "%20".join(species_name.split())
    protein = "%20".join(prot_name.split())
    url = "https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)"\
          "%20AND%20(organism_name:%22{}%22)%20AND%20(protein_name:%22{}%22)"\
          "%20AND%20(structure_3d:true)".format(species, protein)
    with urllib.request.urlopen(url) as hstream:
        data = json.load(hstream)
    if len(data["results"]) == 0:
        logger.debug("The search for a candidate with species {} and name {} did "
                     "not produce any Reviewed entry.".format(species_name, prot_name))
        return None
    if len(data["results"]) > 1:
        logger.debug("The search for a candidate with species {} and name {} "
                     "produced more than one result.".format(species_name, prot_name))
        return None
    logger.debug("Found ac for {} from {} => {}".format(prot_name, species_name,
                 data["results"][0]["primaryAccession"]))
    return data["results"][0]["primaryAccession"]


def search_by_sp_and_name_BCK(species_name, prot_name):
    species = "+".join(species_name.split())
    protein = "+".join(prot_name.split())
    url = 'https://www.uniprot.org/uniprot/?query=organism%3A%22' + species + \
          '%22+name%3A%22' + protein + '%22+database%3A(type%3Apdb)&sort=score'
    document = ""
    tried = 0
    while tried < 3:
        try:
            with urllib.request.urlopen(url) as hstream:
                stream = hstream.read()
                stream = stream.decode("utf-8")
                for line in stream:
                    document += line
                break
        except http.client.RemoteDisconnected:
            tried += 1
            logger.error("ERROR in http connection when searching {} from {}."
                         .format(species_name, prot_name))
        except urllib.error.HTTPError:
            tried += 1
            logger.error("ERROR in http connection when searching {} from {}."
                         .format(species_name, prot_name))
    soup = BeautifulSoup(document, "lxml")
    results_table = soup.find_all("table")
    table = None
    for result in results_table:
        if result.attrs["class"][0] == "grid" and result.attrs["id"] == "results":
            table = result
    if table is not None:
        tbody = result.find("tbody")
        rows_results = tbody.find_all("tr")
        if len(rows_results) == 1:
            for row in rows_results:
                tds_in_row = row.find_all("td")
                for td in tds_in_row:
                    div = td.find("div")
                    if div is not None:
                        try:
                            if div.attrs["title"] == "Reviewed (Swiss-Prot)":
                                logger.error("Found ac for %s from %s => %s"
                                             % (prot_name, species_name,
                                                     row.attrs["id"]))
                                return row.attrs["id"]
                            else:
                                logger.error("The search for a candidate with species %s and name %s did "
                                             "not produce any Reviewed entry." % (species_name, prot_name))
                                return None
                        except KeyError:
                            logger.error("An unexpected error occurred handling the protein %s from "
                                         "%s. Please check the page %s manually." % (prot_name, species_name, url))
                            return None
        else:
            logger.error("The search for a candidate with species %s and name "
                         "%s produced more than one result." % (species_name, prot_name))
            return None
    else:
        logger.error("The search for a candidate with species %s and name %s did "
                     "not produce any results." % (species_name, prot_name))
        return None

