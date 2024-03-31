#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
See online documentation at:
    https://openbabel.org/docs/dev/Command-line_tools/babel.html
"""

import os

from external import process
from libs.configuration import ENV


def convert(in_file, out_file, in_format=None, out_format=None):
    args = []
    if in_format is not None:
        args.extend(["-i", in_format])
    args.append(in_file)
    if out_format is not None:
        args.extend(["-o", out_format])
    args.extend(["-O", out_file])

    _exec_babel(args)

    if not os.path.exists(out_file):
        raise RuntimeError("Failed to convert molecule with OpenBabel.")


def _exec_babel(args):
    process.execute([ENV["babel"]] + args, shell=ENV["babel_shell"])
