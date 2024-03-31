#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from external import process
from libs.configuration import ENV
from libs import paths


def execute(input_path, file_name):
    if ENV["clustal_omega"] is None:
        output_path = os.path.join(paths.manual(), file_name)
        _manual_input(input_path, output_path)
        return output_path
    else:
        output_path = os.path.join(paths.temp(), file_name)
        _execute_local(input_path, output_path)
        return output_path


def _manual_input(input_path, output_path):
    if not os.path.exists(output_path):
        message = "Missing Clustalo result file. \nPlease take input from " \
                  "'{}' and use it with Clustalo service.\n" \
                  "https://toolkit.tuebingen.mpg.de/#/tools/clustalo \n" \
                  "Save the output as '{}'".format(input_path, output_path)
        raise RuntimeError(message)


def _execute_local(input_path, output_path):
    """
    output overwrite (--force) used to enable iterative running of 00_prepare_knowledge_base.py
    otherwise we will have an runtimeerror if the output exists
    """
    command = [ENV["clustal_omega"], "-i", input_path, "-o", output_path, "--force"]
    process.execute(command, shell=ENV["clustal_omega_shell"])
