#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from libs.data_model import SubstrateEntry

_ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")

_DATA_DIR = None


def as_absolute(path):
    """
    Resolve given relative path using project root.
    """
    return os.path.join(_ROOT, path)


def set_task_name(task_id):
    global _DATA_DIR
    _DATA_DIR = os.path.join(_ROOT, "data", task_id)
    os.makedirs(_DATA_DIR, exist_ok=True)


def task_root():
    return _DATA_DIR


def logs():
    path = os.path.join(_DATA_DIR, "logs")
    os.makedirs(path, exist_ok=True)
    return path


def manual():
    path = os.path.join(_DATA_DIR, "manual")
    os.makedirs(path, exist_ok=True)
    return path


def filters():
    path = os.path.join(_DATA_DIR, "filters")
    os.makedirs(path, exist_ok=True)
    return path


def substrates_clusters():
    path = os.path.join(_DATA_DIR, "clusters")
    os.makedirs(path, exist_ok=True)
    return path


def candidates():
    path = os.path.join(_DATA_DIR, "candidates")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_sdf():
    path = os.path.join(candidates(), "sdf")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_pdbqt():
    path = os.path.join(candidates(), "pdbqt")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_docking():
    path = os.path.join(candidates(), "docking")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_pubchem():
    path = os.path.join(candidates(), "pubchem")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_info():
    return os.path.join(candidates(), "info.json")


def candidates_filtered():
    path = os.path.join(_DATA_DIR, "candidates_filtered")
    os.makedirs(path, exist_ok=True)
    return path


def candidates_vina_template():
    return os.path.join(candidates(), "vina_template.txt")


def target_protein():
    return os.path.join(candidates(), "protein.pdbqt")


def substrates_eval():
    return os.path.join(_DATA_DIR, "substrates-evaluation.jsonl")


def get_candidates_eval_id(eval_filename):
    if "candidates-evaluation." in eval_filename and eval_filename.endswith(".jsonl.gz"):
        return os.path.basename(eval_filename)[len("candidates-evaluation."):-len(".jsonl.gz")]


def candidates_eval(file_id=None):
    if file_id is not None:
        return os.path.join(_DATA_DIR, "candidates-evaluation." + str(file_id) + ".jsonl.gz")
    else:
        file_list = []
        for eval_filename in os.listdir(_DATA_DIR):
            if "candidates-evaluation." in eval_filename and eval_filename.endswith(".jsonl.gz"):
                file_list.append(os.path.join(_DATA_DIR, eval_filename))
        return file_list


def knowledge_base_path():
    return os.path.join(_DATA_DIR, "knowledge_base")


def substrate_sdf_path(substrate: SubstrateEntry):
    return os.path.join(_DATA_DIR, "substrates",
                        file_name_for_substrate(substrate) + ".sdf")


def file_name_for_substrate(substrate):
    database = substrate.database().first().value()
    if database == "pubchem":
        cid = substrate.prop("cid").first().value()
        return "pubchem_{}".format(str(cid))
    else:
        raise RuntimeError("Unknown database", database)


def substrates():
    path = os.path.join(_DATA_DIR, "substrates")
    os.makedirs(path, exist_ok=True)
    return path


def substrate_pdbqt_path(substrate: SubstrateEntry):
    return os.path.join(_DATA_DIR, "substrates",
                        file_name_for_substrate(substrate) + ".pdbqt")


def protein_pdbqt_path(mapping, mcsa):
    return os.path.join(_DATA_DIR, "proteins",
                        file_name_for_mapping(mapping, mcsa) + ".pdbqt")


def file_name_for_mapping(mapping, mcsa):
    source_pdb = mapping.source_pdb().first().value()
    source_chain = mapping.source_chain().first().value()
    target_pdb = mapping.target_pdb().first().value()
    target_chain = mapping.target_chain().first().value()
    return "{}-{}_{}-{}_{}".format(
        source_pdb, source_chain, target_pdb, target_chain, mcsa.mcsa_id())


def database_directory():
    return os.path.join(_DATA_DIR, "knowledge_base")


def temp():
    path = os.path.join(_DATA_DIR, "temp")
    os.makedirs(path, exist_ok=True)
    return path


def proteins():
    path = os.path.join(_DATA_DIR, "proteins")
    os.makedirs(path, exist_ok=True)
    return path


def complexes():
    path = os.path.join(_DATA_DIR, "complexes")
    os.makedirs(path, exist_ok=True)
    return path


def complexes_pdbqt_path(name):
    return os.path.join(_DATA_DIR, "complexes", name + ".pdbqt")


def accessibility_visualisation(fragment):
    file_name = "{}-{}-{}".format(
        str(fragment["name"]),
        str(fragment["model_index"]),
        str(fragment["fragment_index"]))
    return os.path.join(_DATA_DIR, "visualisation", file_name)


def task_info():
    return os.path.join(_DATA_DIR, "task.json")


def get_subdir_filename_from_filepath(file_path):
    tmp_path, filename = os.path.split(file_path)
    subdir = os.path.basename(tmp_path)
    return subdir, filename
