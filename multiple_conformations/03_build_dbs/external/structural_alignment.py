#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Align given file (from_file) to match position of template file (to_file).
"""

import os

from external import process
from libs import paths
from libs.utils import remove_existing
from libs.configuration import ENV


def align_pdb_files(from_file, to_file, residue_mapping, output_path):
    script = _create_align_fit_script(
        from_file, to_file, residue_mapping, output_path)

    script_path = os.path.join(paths.temp(), "structural_alignment.py")
    with open(script_path, "w") as output_stream:
        output_stream.write(script)

    process.execute([ENV["python_pymol"], script_path],
                    shell=ENV["python_pymol_shell"])

    if not os.path.exists(to_file):
        raise RuntimeError("Failed to align structure (no output). You can check the script at: " + script_path)

    remove_existing(script_path)


def _create_align_fit_script(from_file, to_file, mapping, output_path):
    script = "\n"
    script += "import pymol\n"
    script += "from pymol import cmd\n"
    script += "pymol.pymol_argv = [\"pymol\", \"-qc\"]\n"
    script += "pymol.finish_launching()\n"
    script += "\n"
    script += "cmd.load(\"" + from_file.replace("\\", "/") + \
              "\", \"f\", 0, \"pdb\")\n"
    script += "cmd.load(\"" + to_file.replace("\\", "/") + \
              "\", \"t\", 0, \"pdb\")\n"
    script += "\n"

    f_select = " or ".join(
        "/f//{}/{}".format(item["from"]["chain"], item["from"]["res"])
        for item in mapping)
    script += "cmd.select(\"f_sel\", \"" + f_select + "\")\n"

    t_select = " or ".join(
        "/t//{}/{}".format(item["to"]["chain"], item["to"]["res"])
        for item in mapping)
    script += "cmd.select(\"t_sel\", \"" + t_select + "\")\n"

    script += "\n"

    script += "cmd.align(\"f_sel\", \"t_sel\")\n"

    script += "\n"
    script += "cmd.select(\"out\", \"/f////\")\n"
    script += "cmd.save(\"" + output_path.replace("\\", "/") + "\", \"out\")\n"
    return script
