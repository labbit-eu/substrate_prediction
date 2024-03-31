#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from libs import paths
import logging
import subprocess

logger = logging.getLogger(__name__)


def execute(command, shell, cwd=None):
    logger.debug("Executing command %s", " ".join(map(str, command)))
    logger.debug("  in shell: %s with cwd: '%s'", str(shell), str(cwd))
    _check_arguments(shell)
    _execute_command(command, shell, cwd)


def _check_arguments(shell):
    if not isinstance(shell, bool):
        raise RuntimeError("Parameter 'shell' must be boolean.")


def _execute_command(command, shell, cwd):
    with subprocess.Popen(
            command,
            cwd=cwd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE) as process:
        out, err = process.communicate()
        errcode = process.returncode
        if "clustalw2" in command[0]:
            logfile = os.path.join(paths.temp(), "similarity.log")
            with open(logfile, "w") as outstream:
                outstream.write(out.decode("utf-8"))
        elif "osra" in command[0]:
            logfile = os.path.join(paths.temp(), "smiles.log")
            with open(logfile, "w") as outstream:
                outstream.write(out.decode("utf-8"))
        # elif "molconvert" in command[0]:
        #     logfile = os.path.join(paths.temp(), "compound.log")
        #     with open(logfile, "wb") as outstream:
        #         outstream.write(out)
        #         outstream.write(err)
    if errcode == 0:
        return
    message = "Shell command failed: \n"
    message += " ".join(command)
    message += "\nerror code:"
    message += str(errcode)
    if out is not None:
        message += "\nstdout:\n"
        message += out.decode("utf-8")
    if err is not None:
        message += "\nstderr:\n"
        message += err.decode("utf-8")
    raise RuntimeError(message)
