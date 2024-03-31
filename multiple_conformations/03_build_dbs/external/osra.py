#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2018 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#

from external import process
from libs.configuration import ENV


def execute(input_path):
    command = [ENV["osra"], "-f smi", input_path]
    process.execute(command, shell=ENV["osra_shell"])


