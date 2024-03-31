#!/usr/bin/env python
# -*- coding: utf-8 -*-

from libs import utils


def download_sdf_as_path(chebi):
    url = "https://www.ebi.ac.uk/chebi/saveStructure.do?" \
          "sdf=true&chebiId=" + str(chebi)
    return utils.fetch_url(chebi, url, "chebi-sdf")
