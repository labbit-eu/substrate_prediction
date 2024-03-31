#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from external import process
from libs.configuration import ENV


def prepare_ligand(in_path, out_path):
    """
    http://autodock.scripps.edu/faqs-help/how-to/how-to-prepare-a-ligand-file-for-autodock4
    """
    script = "Utilities24/prepare_ligand4.py"

    _run_auto_dock_tools_script(script, ["-l", in_path, "-o", out_path])

    if not os.path.exists(out_path):
        raise RuntimeError("Failed to prepare ligand with MglTools.")


def _run_auto_dock_tools_script(script, arguments=None):
    script_path = os.path.join(ENV["mgl"], script)
    command = [ENV["python_mgl"], script_path] + (arguments or [])
    process.execute(command, ENV["python_mgl_shell"])


def prepare_receptor(in_path, out_path):
    """
    http://autodock.scripps.edu/faqs-help/how-to/how-to-prepare-a-receptor-file-for-autodock4
    """
    script = "Utilities24/prepare_receptor4.py"
    # arguments minimizing the thorough cleaning that is provided by external/standardize_protein_pdb_file
    args = ["-A", "checkhydrogens",
            "-U", "nphs_lps"]
    
    _run_auto_dock_tools_script(script, ["-r", in_path, "-o", out_path] + args)

    if not os.path.exists(out_path):
        raise RuntimeError("Failed to prepare protein with MglTools.")
