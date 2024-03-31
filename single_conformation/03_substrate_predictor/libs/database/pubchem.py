#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilize pubchempy.py as a third party library.
"""

import json
import shutil
import logging

from libs.database import pubchempy
from libs import utils
from libs.utils import ObjectCache

logger = logging.getLogger(__name__)


class PubchemCompound(object):

    def __init__(self, cid):
        self.cid = cid
        self.smiles = None
        self.name = None


class PubchemSearchByNameCompoundFactory(object):
    cache_name = "pubchem-compound"

    def create(self, query):
        records = []
        for compound in pubchempy.get_compounds(query, "name"):
            record = PubchemCompound(compound.cid)
            record.smiles = compound.isomeric_smiles
            record.name = compound.iupac_name
            records.append(record)
        return records

    @staticmethod
    def to_json(data):
        return [item.__dict__ for item in data]

    @staticmethod
    def from_json(data):
        records = []
        for item in data:
            record = PubchemCompound(item["cid"])
            record.smiles = item["smiles"]
            record.name = item["name"]
            records.append(record)
        return records


def query_by_name(name):
    pubchem_search_factory = PubchemSearchByNameCompoundFactory()
    with ObjectCache(pubchem_search_factory, name) as molecules:
        output = molecules
    return output


def download_3d_sdf(cid, path):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + \
          str(cid) + "/record/SDF/?record_type=3d"
    cache_path = utils.fetch_url(cid, url, "pubchem-compound-3d-sdf")
    shutil.copy(cache_path, path)


def download_3d_sdf_as_path(cid):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + \
          str(cid) + "/record/SDF/?record_type=3d"
    return utils.fetch_url(cid, url, "pubchem-compound-3d-sdf")


def download_similar_3d_sdf_as_path(cid, threshold):
    output = []
    for similar_cid in fetch_similar_cid(cid, threshold):
        try:
            output.append(download_3d_sdf_as_path(similar_cid))
        except Exception:
            """
            There might be compounds without structure. An example is 88372870, 
            that utilize structure of 10898. The first one mentioned have 
            extra free As atom. 
            """
            logger.exception(
                "Failed to download structure for: %s", similar_cid)
    return output


def fetch_similar_cid(cid, threshold):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/" \
          "fastsimilarity_2d/cid/{}/cids/JSON?Threshold={}" \
        .format(cid, threshold)
    id_ = "{}-{}".format(cid, threshold)
    cache_name = "pubchem-fastsimilarity-2d"
    with open(utils.fetch_url(id_, url, cache_name)) as input_stream:
        content = json.load(input_stream)
    return content["IdentifierList"]["CID"]
