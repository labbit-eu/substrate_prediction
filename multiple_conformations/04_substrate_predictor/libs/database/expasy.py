#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Data source URL:
    https://enzyme.expasy.org/EC/{ec number}.txt
"""

from libs.utils import UrlCache


class ExpasyEntry(object):

    def __init__(self, ec_number):
        self.ec_number = ec_number
        self.uniprot_acs = []

    def __str__(self):
        return self.ec_number


def query_entry(ec_number):
    url = "https://enzyme.expasy.org/EC/" + ec_number + ".txt"
    entry = ExpasyEntry(ec_number)
    with UrlCache(ec_number, url, "expasy") as stream:
        for line in stream:
            if line.startswith("DR"):
                entry.uniprot_acs.extend(_parse_dr_line(line))
    return entry


def _parse_dr_line(line):
    unipro_acs = []
    for chunk in line[3:].split(";"):
        chunk = [x.strip() for x in chunk.split(",")]
        if not len(chunk) == 2:
            continue
        unipro_acs.append(chunk[0])
    return unipro_acs
