#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import argparse
import logging
import math
import configparser

logger = logging.getLogger(__name__)

ENV = {
    "vina": None,
    "vina_shell": False,
    "python_mgl": None,
    "python_mgl_shell": False,
    "mgl": None,
    "babel": None,
    "babel_shell": False,
    "python_pymol": None,
    "python_pymol_shell": False,
    "clustal_omega": None,
    "clustal_omega_shell": False,
    # Path to BRENDA database dump downloaded from:
    # http://www.brenda-enzymes.org/download_brenda_without_registration.php
    "console_log_level": "INFO",
    "max_candidates_in_dir": 5000,
    "num_parallel_cpu": None,
    "cache_directory": "./cache/",
    "pubchem": {
                  # Maximum should be 499999, but lower values make it less painful
                  # if something goes wrong as we lost less data.
                  "batch_size": 200000
              }
}

REMOTE = {
    "sleep_time": 300,
    "selfsustained": True,
    "update_template_file": False,
    "vina_dir": None,
    "batch_size": None,
    "batch_number": None,
    "docking_dataname": "run_docking",
    "host": {
        "hostname": None,
        "username": None,
        "port": 22,
        "submit_dir": None,
        "ssh_key_path": None
    },
    "queue": {
        "manager": None,
        "name": None,
        "mem": 500,
        "tmpdir": None,
        "shell_path": None,
        "template_file": None
    }
}


# _ATOM_RADIUS values are comming from AutoDock vina
_ATOM_RADIUS = {
    "C": 2.00000,
    "A": 2.00000,
    "N": 1.75000,
    "O": 1.60000,
    "P": 2.10000,
    "S": 2.00000,
    "H": 1.00000,
    "F": 1.54500,
    "I": 2.36000,
    "NA": 1.75000,
    "OA": 1.60000,
    "SA": 2.00000,
    "HD": 1.00000,
    "Mg": 0.65000,
    "Mn": 0.65000,
    "Zn": 0.74000,
    "Ca": 0.99000,
    "Fe": 0.65000,
    "Cl": 2.04500,
    "Br": 2.16500
}


_SPECIAL_RESIDUES = {
    "ARG": {
        "atoms": ["NH1", "NH2"],
        "atoms2": ["1HH1", "2HH1"],
        "atoms3": ["1HH2", "2HH2"]
    },
    "ASN": {
        "HD21": "HD2",
        "HD22": "HD2",
        "1HD2": "HD2",
        "2HD2": "HD2"
    },
    "ASP": {
        "OD1": "OD",
        "OD2": "OD"
    },
    "GLN": {
        "atoms": ["1HE2", "2HE2"]
    },
    "GLU": {
        "OE1": "OE",
        "OE2": "OE"
    },
    "LEU": {
        "atoms": ["CD1", "CD2"]
    },
    "LYS": {
        "atoms": ["HZ1", "HZ2", "HZ3"],
        "HZ1": "HZ",
        "HZ2": "HZ",
        "HZ3": "HZ"
    },
    "VAL": {
        "atoms": ["CG1", "CG2"]
    }
}


# here we now follow the standard PDB residue naming convention: ftp://ftp.wwpdb.org/pub/pdb/data/monomers
_RESIDUE_TYPE_FUNCTION = {
    ("ASN", "electrostatic stabiliser"): {
        "atoms": ["HD21", "HD22"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("ASP", "covalent catalysis"): {
        "atoms": ["OD1", "OD2"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("ASP", "electrostatic stabiliser"): {
        "atoms": ["OD1", "OD2"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("CYS", "electrostatic stabiliser"): {
        "atoms": ["SG"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("CYS", "metal ligand"): {
        "atoms": ["SG"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("GLN", "electrostatic stabiliser"): {
        "atoms": ["OE1", "NE2"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("GLU", "electrostatic stabiliser"): {
        "atoms": ["OE1", "OE2"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("HIS", "proton shuttle (general acid/base)"): {
        "atoms": ["ND1", "NE2"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("HIS", "metal ligand"): {
        "atoms": ["ND1", "NE2"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("HIS", "electrostatic stabiliser"): {
        "atoms": ["ND1", "NE2", "HD1", "HE2"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("LEU", "electrostatic stabiliser"): {
        "atoms": ["N"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("LYS", "electrostatic stabiliser"): {
        "atoms": ["NZ", "HZ1", "HZ2", "HZ3"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("LYS", "modifies pKa"): {
        "atoms": ["NZ"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("PHE", "steric role"): {
        "atoms": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("PRO", "steric role"): {
        "atoms": ["CB", "CG", "CD"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("SER", "covalent catalysis"): {
        "atoms": ["OG"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("SER", "proton shuttle (general acid/base)"): {
        "atoms": ["OG"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("THR", "electrostatic stabiliser"): {
        "atoms": ["OG1"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("THR", "steric role"): {
        "atoms": ["CB", "OG1", "CG2"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("THR", "proton shuttle (general acid/base)"): {
        "atoms": ["CB", "OG1", "CG2"],
        "threshold": 5,
        "tolerance": 0.05,
        "accessibility_threshold": 0.5
    },
    ("TRP", "electrostatic stabiliser"): {
        "atoms": ["HE1"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("TYR", "proton donor"): {
        "atoms": ["OH"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("TYR", "proton acceptor"): {
        "atoms": ["OH"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("TYR", "electrostatic stabiliser"): {
        "atoms": ["HH"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("NDP", "proton donor"): {
        "atoms": ["C4N"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("NDP", "electrostatic stabiliser"): {
        "atoms": ["C4N"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("NAP", "electrostatic stabiliser"): {
        "atoms": ["C4N"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    },
    ("NPH", "electrostatic stabiliser"): {
        "atoms": ["C4N"],
        "threshold": 6,
        "tolerance": 0.5,
        "accessibility_threshold": 0.5
    }
}

"""
Used as a threshold to evaluate whether a general type of interaction 
can take a place between atoms.

Used in docking_evaluation to check, whether atoms of opposite charge
are close enough.
"""
_ACTIVITY_THRESHOLD = 8

_ALLOWED_ATOMS = [
    "C", "N", "O", "P", "S", "H", "F", "I", "Mg", "Mn", "Zn", "Ca", "Fe",
    "Cl", "Br"
]

_ACCESSIBILITY = {
        "sphere_definition": [
            # [vertical angle, horizontal angle [start, end, step], weight]
            [math.radians(0), [0, math.pi, math.pi], 0.4],
            [math.radians(15), [0, 2 * math.pi, math.radians(60)], 0.3],
            [math.radians(30), [0, 2 * math.pi, math.radians(45)], 0.2],
            [math.radians(45), [0, 2 * math.pi, math.radians(30)], 0.1]
        ],
        "pre_calculated_sphere": [
            [[ 0.000, 1.000, 0.000], 0.4],
            [[ 0.259, 0.966, 0.000], 0.05],
            [[ 0.129, 0.966, 0.224], 0.05],
            [[-0.129, 0.966, 0.224], 0.05],
            [[-0.259, 0.966, 0.000], 0.05],
            [[-0.129, 0.966,-0.224], 0.05],
            [[ 0.129, 0.966,-0.224], 0.05],
            [[ 0.500, 0.866, 0.000], 0.025],
            [[ 0.354, 0.866, 0.354], 0.025],
            [[ 0.000, 0.866, 0.500], 0.025],
            [[-0.354, 0.866, 0.354], 0.025],
            [[-0.500, 0.866, 0.000], 0.025],
            [[-0.354, 0.866,-0.354], 0.025],
            [[-0.000, 0.866,-0.500], 0.025],
            [[ 0.354, 0.866,-0.354], 0.025],
            [[ 0.707, 0.707, 0.000], 0.0083],
            [[ 0.612, 0.707, 0.354], 0.0083],
            [[ 0.354, 0.707, 0.612], 0.0083],
            [[ 0.000, 0.707, 0.707], 0.0083],
            [[-0.354, 0.707, 0.612], 0.0083],
            [[-0.612, 0.707, 0.354], 0.0083],
            [[-0.707, 0.707, 0.000], 0.0083],
            [[-0.612, 0.707,-0.354], 0.0083],
            [[-0.354, 0.707,-0.612], 0.0083],
            [[-0.000, 0.707,-0.707], 0.0083],
            [[ 0.354, 0.707,-0.612], 0.0083],
            [[ 0.612, 0.707,-0.354], 0.0083]
        ]
    }

CHEM = {
    "atom_radius": _ATOM_RADIUS,
    "residues": _RESIDUE_TYPE_FUNCTION,
    "activity_threshold": _ACTIVITY_THRESHOLD,
    "allowed_atoms": _ALLOWED_ATOMS,
    "accessibility": _ACCESSIBILITY
}

TASK = {
    "name": None,
    "knowledge_base": {
        "uniprot_search_remote": True,
        "uniprot_search_evalue": 0.001,
        "ec_number": None,
        "ec_support_threshold": 0.8,
        "pdb": None,
        "protein_chain": None,
        "hetero_residues_to_keep": None
    },
    "docking": {
        "box": None,
        "cpu": 4,
        "exhaustiveness": 8,
        "seed": 849,
        "num_modes": 100,
        "energy_range": 10
    },
    "substrates": {
        "cluster": {
            "merge_threshold": 1.0,
            "method": "complete",
            "min_cluster_size": 3
        },
        "cluster_filter": {
            "input": "^[0-9]+$",
            "backup": None
        },
        "accessibility": {
            "input": "^[0-9]+$",
            "backup": None,
            "save_visualisation": True
        },
        "filter_create": {
            "input": "^[0-9]+$",
            "relax_threshold": 0.075,
            "distance_type": "single",
            "accessibility_type": "single",
            "normalization": "range"
        },
        "group_filters": {
            "input": "^[0-9]+$",
            "backup": None
        }
    },
    "candidates": {
        "docking_mode": "internal",
        "prepare": {
            "max_molecules": 100000,
            "logP": [-10, 8],
            "MW": [0, 500],
            "TORSDOF": [0, 15],
            "manual": None
        },
        "filter_distances": {
            "input": "^[0-9]+$"
        },
        "accessibility": {
            "input": "^[0-9]+$",
            "backup": None,
            "save_visualisation": False
        },
        "filter_accessibility": {
            "input": "^[0-9]+$",
            "backup": None
        },
        "print_detail": {
            "input": "^[0-9]+$"
        },
        "print_summary": {
            "input": "^[0-9]+$"
        },
        "cluster_plot_info": {
            "input": "^[0-9]+$"
        },
        "letters": ["_A"]
    },
    "distances": {
    "selected": ("resid", "role", "atname", "atname")
    }
}


def initialize_configuration(module_file=None):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    args, unknown_args = _read_args()

    for path in args.config:
        config_path = os.path.join(this_dir, "..", path)
        _load_configuration(config_path)
        
    if module_file is not None:
        for key_name in TASK.keys():
            if key_name in os.path.basename(module_file):
                _merge_args_to_configuration(key_name, unknown_args)

    if TASK["knowledge_base"]["pdb"] is None:
        raise RuntimeError("PDB-ID of target protein was not specified in the parameter \"knowledge_base.pdb\" in "
                           "[TASK] section of the configuration file.")

    # TASK["knowledge_base"]["pdb"] = TASK["knowledge_base"]["pdb"].upper()

    ENV["console_log_level"] = _logging_level_from_string(ENV["console_log_level"])
    ENV["num_parallel_cpu"] = get_num_cpu(ENV["num_parallel_cpu"])
    
    if TASK["candidates"]["docking_mode"] == "remote":
        _test_remote_configuration()

    if isinstance(TASK["candidates"]["letters"], str):
        TASK["candidates"]["letters"] = [l.strip() for l in TASK["candidates"]["letters"].split(",")]

    return args


def _test_remote_configuration():
    REMOTE["queue"]["shell_path"] = _get_batch_shell_path(REMOTE["queue"]["shell_path"])
    if "bash" not in REMOTE["queue"]["shell_path"]:
        raise RuntimeError("Selected shell ({}) is not supported.  Please alter the parameter \"queue.hell_path\" in "
                           "[REMOTE] section of the configuration file.".format(REMOTE["queue"]["shell_path"]))

    if REMOTE["queue"]["manager"] is None:
        raise RuntimeError("The type of queue manager for remote submission is not set. Please specify the name in "
                           "the parameter \"queue.manager\" in [REMOTE] section of the configuration file.")
    else:
        REMOTE["queue"]["manager"] = REMOTE["queue"]["manager"].lower()
        if REMOTE["queue"]["manager"] not in ["slurm", "torque"]:
            raise RuntimeError("Selected queue manager ({}) is not supported. Please alter the parameter \"queue."
                               "manager\" in [REMOTE] section of "
                               "the configuration file.".format(REMOTE["queue"]["manager"]))

    parameters = [("queue", "name", "The name of queue for remote submission is not set."),
                  ("queue", "tmpdir", "The temporary directory used during remote execution of jobs is not set."),
                  ("host", "hostname", "The hostname for remote submission is not set."),
                  ("host", "username", "The username for remote submission is not set."),
                  ("host", "submit_dir", "The root directory for remote submission on host is not set."),
                  ("host", "ssh_key_path", "The path to ssh_key enabling access to remote host is not set.")]

    for param in parameters:
        if REMOTE[param[0]][param[1]] is None:
            raise RuntimeError("{} Please specify it in the parameter \"{}.{}\" in [REMOTE] section of the "
                               "configuration file.".format(param[2], param[0], param[1]))

    if (REMOTE["batch_size"] is None) and (REMOTE["batch_number"] is None):
        message = "For remote candidate screening, one of the parameters \"batch_size\" or " \
                  "\"batch_number\" have to be specified in [REMOTE] section of the configuration file."
        raise RuntimeError(message)

    if (REMOTE["batch_size"] is not None) and (REMOTE["batch_number"] is not None):
        message = "For remote candidate screening, only one of the parameters \"batch_size\" and " \
                  "\"batch_number\" has to be specified in [REMOTE] section of the configuration file."
        raise RuntimeError(message)


def _read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", dest="config",
                        action="append", help="Configuration path.",
                        required=True)
    parser.add_argument("-f", "--filters", dest="filters",
                        action="store", help="Selected filters to test against.",
                        required=False)
    return parser.parse_known_args()


def _merge_args_to_configuration(module_name, args):
    for index in range(0, len(args)):
        if not args[index].startswith("--"):
            continue
        key = args[index][2:]
        if index + 1 >= len(args):
            raise RuntimeError("Missing argument for: " + key)

        value = args[index + 1]
        TASK[module_name][key] = value


def _load_configuration(path):
    paths = list(_resolve_configuration_path(path))
    if len(paths) == 0:
        raise RuntimeError("Missing configuration file: " + path)

    for path in paths:
        logger.info("Loading configuration from: %s", path)
        if path.lower().endswith(".json"):
            _load_configuration_json(path)
        elif path.lower().endswith(".ini"):
            _load_configuration_ini(path)
        else:
            raise RuntimeError("Invalid configuration file type: " + path)


def _resolve_configuration_path(path):
    if os.path.exists(path):
        yield path
    if os.path.exists(path + ".json"):
        yield path + ".json"
    if os.path.exists(path + ".ini"):
        yield path + ".ini"


def _load_configuration_json(path):
    with open(path) as in_stream:
        config = json.load(in_stream)
    try:
        _merge(config, {"task": TASK, "env": ENV, "remote": REMOTE, "chem": CHEM})
    except RuntimeError:
        print("Invalid configuration: " + path)
        raise


def _load_configuration_ini(path):
    float_array_paths = {"docking.box",
                         "candidates.prepare.logP",
                         "candidates.prepare.MW",
                         "candidates.prepare.TORSDOF"
                         }
    str_array_paths = {"knowledge_base.hetero_residues_to_keep"}
    boolean_paths = {"vina_shell",
                     "python_mgl_shell",
                     "babel_shell",
                     "knowledge_base.uniprot_search_remote",
                     "python_pymol_shell",
                     "clustal_omega_shell",
                     "update_template_file",
                     "selfsustained",
                     "knowledge_base.uniprot_search_remote",
                     "substrates.accessibility.save_visualisation",
                     "candidates.accessibility.save_visualisation"
                     }
    integer_paths = {"max_candidates_in_dir",
                     "num_parallel_cpu",
                     "pubchem.batch_size",
                     "sleep_time",
                     "batch_size",
                     "batch_number"
                     "host.port",
                     "queue.mem",
                     "docking.cpu",
                     "docking.exhaustiveness",
                     "docking.seed",
                     "docking.num_modes",
                     "docking.energy_range",
                     "substrates.cluster.min_cluster_size",
                     "candidates.prepare.max_molecules"
                     }
    float_paths = {"knowledge_base.ec_support_threshold",
                   "knowledge_base.uniprot_search_evalue",
                   "substrates.cluster.merge_threshold",
                   "filter_create.relax_threshold"
                   }
    file_paths = {"vina",
                  "python_mgl",
                  "mgl",
                  "babel",
                  "python_pymol",
                  "clustal_omega",
                  "queue.template_file"
                  }
    others = {"distances.selected",
              "distances.selected_residues"
              }

    config = configparser.ConfigParser()
    config.read(path)
    config_sections = {"task": TASK, "env": ENV, "remote": REMOTE, "chem": CHEM}
    for section in config.sections():
        root = config_sections[section]
        for name, value in config[section].items():
            # Get object to set value to.
            path = name.split(".")
            dictionary = root
            for key in path[:-1]:
                dictionary = dictionary[key]
            # Prepare value.
            if name in float_array_paths:
                value = [_test_float(name, val, section) for val in value.split(",")]
            elif name in boolean_paths:
                value = _test_boolean(name, value, section)
            elif name in integer_paths:
                value = _test_int(name, value, section)
            elif name in float_paths:
                value = _test_float(name, value, section)
            elif name in str_array_paths:
                value = [str(val) for val in value.split(",")]
            elif name in file_paths:
                _test_file_paths(name, value, section)
            elif name in others:
                value = _test_others(name, value, section)
            dictionary[path[-1]] = value


def _test_others(name, value, section):
    keys = []
    if name == "distances.selected":
        try:
            for line in value.split("\n"):
                chunks = [s.strip() for s in line.split(",")]
                assert len(chunks) == 4
                keys.append((chunks[0], chunks[1], chunks[2], chunks[3]))
        except AssertionError:
            raise AssertionError("Key for distance {} is incomplete".format(line))
    elif name == "distances.selected_residues":
        for line in value.split("\n"):
            chunks = [s.strip() for s in line.split()]
            keys.extend(chunks)
    return keys


def _test_int(name, value, section):
    try:
        return int(value)
    except ValueError:
        if float(value).is_integer():
            return int(float(value))
        else:
            raise ValueError("Invalid value \"{}\" of parameter \"{}\" in the [{}] section of the configuration file. "
                             "It has to be an integer.".format(value, name, section.upper()))


def _test_float(name, value, section):
    try:
        return float(value)
    except ValueError:
        raise ValueError("Invalid value \"{}\" of parameter \"{}\" in the [{}] section of the configuration file. "
                         "It has to be float.".format(value, name, section.upper()))


def _test_file_paths(parameter, file_path, section):
    file_name = os.path.basename(file_path)
    dir_name = os.path.dirname(file_path)
    if not os.path.exists(file_path):
        raise RuntimeError("Could not find {} file in following location: {}. Please verify that you have properly "
                           "set the path in the parameter \"{}\" in [{}] section of "
                           "the configuration file".format(file_name, dir_name, parameter, section.upper()))


def _test_boolean(name, value, section):
    lower_value = value.lower()
    if lower_value == "true" or lower_value == "1":
        return True
    if lower_value == "false" or lower_value == "0":
        return False
    raise ValueError("Invalid value \"{}\" of parameter \"{}\" in the [{}] section of the configuration file. It has "
                     "to be boolean (case insensitive false, true or 1, 0).".format(value, name, section.upper()))


def _merge(source, destination):
    for key, value in source.items():
        if key not in destination:
            raise RuntimeError("Unexpected key: " + key + "")
        if isinstance(value, dict):
            _merge(source[key], destination[key])
        else:
            destination[key] = source[key]


def _logging_level_from_string(value):
    value = value.upper()

    mapping = {
        'CRITICAL': logging.CRITICAL,
        'ERROR': logging.ERROR,
        'WARN': logging.WARNING,
        'WARNING': logging.WARNING,
        'INFO': logging.INFO,
        'DEBUG': logging.DEBUG,
        'NOTSET': logging.NOTSET,
    }

    return mapping.get(value, logging.NOTSET)


def get_num_cpu(value):
    cpu_to_use = 1
    cpu_os = os.cpu_count()
    if value is None:
        if cpu_os > 4:
            # many cpus => keep 2 for system
            cpu_to_use = cpu_os - 2
            
        elif cpu_os > 2:
            cpu_to_use = cpu_os - 1
            # 2-4 cpus => keep 1 for system
    else:
        cpu_to_use = value
        if cpu_to_use > cpu_os:
            logger.warning("The number of CPU(%d) specified for parallel processing is larger than the number of "
                           "CPU(%d) reported as available on this computer.\n Consider decreasing this number defined "
                           "in the parameter \"num_parallel_cpu\" in [ENV] section of the configuration file.",
                           cpu_to_use, cpu_os)
        
    if cpu_to_use > 10:
        logger.warning("The number of CPU(%d) specified for parallel processing is quite large, which might actually"
                       " hamper the performance due to IO restrictions.\n Consider decreasing this number defined in "
                       "the parameter \"num_parallel_cpu\" in [ENV] section of the configuration file.", cpu_to_use)

    return cpu_to_use
    
    
def _get_batch_shell_path(value):
    if value is None:
        default_shell = os.environ["SHELL"]
        if default_shell is not None:
            logger.info("The path to shell executable was not defined, using the default shell: %s",
                        os.environ["SHELL"])
            return default_shell
        else:
            raise RuntimeError("The path to shell executable was not defined, and the default shell cannot "
                               "be determined. Please specify the path in the parameter \"queue.shell_path\" "
                               "in [REMOTE] section of the configuration file.")
    else:
        return value    
