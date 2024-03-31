#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2018 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#

from external import process
from libs.configuration import ENV


def execute(smiles, out_path):
    command = [ENV["molconvert"], '"png:w500" -s', '"%s"' % smiles, '-o "%s"' % out_path]
    process.execute(command, shell=True)


