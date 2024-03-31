#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  02_convert_dbs.py
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carse@amu.edu.pl>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import os
import json
import shutil
import logging
import pandas as pd
from libs import utils
from copy import deepcopy
from libs import file_formats
from libs import fragment_extractor
from configparser import ConfigParser


def make_task_json(config):
    os.makedirs(config["OUTPATHS"]["out_path"], exist_ok=True)
    required_mcsa = config["EC_INFO"]["mcsa_id"]
    sdf_path = config["INPATHS"]["fragment_sdf"]
    fragment = fragment_extractor.extract_fragment(sdf_path)
    smiles = file_formats.load_sdf_file_property(sdf_path, "SMILES")
    task_path = config["OUTPATHS"]["path_task"]
    with open(task_path, "w") as out_steam:
        json.dump({
            "mcsa": "@m-csa,{}".format(required_mcsa),
            "fragment": fragment,
            "smiles": smiles
        }, out_steam, indent=2)

def _create_mapping(target_pdb, ref_pdb, mcsa_id, database, uniprot_id, kb_path, chains="A", query=False):
    map_file_name = "@mapping,{}-{}-{}-{}-{}".format(target_pdb, chains, ref_pdb, chains, mcsa_id)
    mapping = {"__init__": [target_pdb, chains, ref_pdb, chains, mcsa_id]}
    mapping["ref_id"] = map_file_name
    mapping["provenance"] = []
    properties = {"target_chain": [{"value": chains, "provenance": []}]}
    properties["offset"] = [{"value": 0, "provenance": []}]
    prop_value = {}
    if query:
        temp_res = [int(res) for res in database["QUERY_PROTEIN"]["active_sites"][target_pdb[:-3]]]
    else:
        temp_res = [int(res) for res in database["UNIPROT"][uniprot_id]["active_sites"][target_pdb]]
    temp_res.sort()
    temp_map = [int(m[1]) for m in database["MAPPING_INFO"]]
    temp_map.sort()
    for i in range(len(temp_map)):
        prop_value[str(temp_res[i])] = temp_map[i]
    properties["fasta_to_mcsa"] = [{"value": prop_value, "provenance": []}]
    properties["is_valid"] = [{"value": True, "provenance": [{"source": "mapping check"}]}]
    properties["m-csa"] = [{"value": "@m-csa,{}".format(mcsa_id), "provenance": []}]
    properties["target_pdb"] = [{"value": ref_pdb, "provenance": []}]
    properties["source_pdb"] = [{"value": target_pdb, "provenance": []}]
    properties["source_chain"] = [{"value": chains, "provenance": []}]
    mapping["properties"] = properties
    with open("{}/mapping/{}.json".format(kb_path, map_file_name), "w") as out_stream:
        json.dump(mapping, out_stream, indent=2)

def _create_protein(pdb_id, ref_pdb, mcsa_id, in_path, out_path):
    in_name = "{}/{}.pdbqt".format(in_path, pdb_id)
    out_name = "{}/{}-A_{}-A_{}.pdbqt".format(out_path, pdb_id, ref_pdb, mcsa_id)
    shutil.copy(in_name, out_name)

def _create_json_complex(pdb_id, subs, uniprot, ref_pdb, mcsa_id, kb_path, chains="A"):
    complex_name = "@complex,{}-{}-@substrate,{}".format(pdb_id, chains, subs)
    map_file_name = "@mapping,{}-{}-{}-{}-{}".format(pdb_id, chains, ref_pdb, chains, mcsa_id)
    kb_complex = {"__init__": [pdb_id, chains, "@substrate,{}".format(subs)]}
    kb_complex["ref_id"] = complex_name
    kb_complex["provenance"] = []
    properties = {"protein": [{"value": "@uniprot,{}".format(uniprot),
                               "provenance": []}],
                  "substrate": [{"value": "@substrate,{}".format(subs),
                                 "provenance": []}],
                  "mapping": [{"value": map_file_name,
                               "provenance": []}]}
    kb_complex["properties"] = properties
    with open("{}/complex/{}.json".format(kb_path, complex_name), "w") as out_stream:
        json.dump(kb_complex, out_stream, indent=2)

def _create_complex(pdb_id, subs, ref_pdb, mcsa_id, in_path, out_path):
    # TODO: now only using one molecular model A for the analysis, change it to
    # TODO: use more chains
    in_name = "{}/{}_@_{}.pdbqt".format(in_path, pdb_id, subs)
    out_name = "{}/{}-A_{}-A_{}_{}.pdbqt".format(out_path, pdb_id, ref_pdb,
                                                 mcsa_id, subs)
    shutil.copy(in_name, out_name)
    cp_name = "{}-A_{}-A_{}_{}".format(pdb_id, ref_pdb, mcsa_id, subs)
    kb_complex = {"name": cp_name}
    ref = {"protein": "@uniprot,NONE",
           "substrate": "@substrate,{}".format(subs),
           "mapping": "@mapping,{}-A-{}-A-{}".format(pdb_id, ref_pdb, mcsa_id)}
    kb_complex["ref"] = ref
    with open("{}/{}.json".format(out_path, cp_name), "w") as out_stream:
        json.dump(kb_complex, out_stream, indent=2)

def _copy_log(pdb_id, subs, ref_pdb, mcsa_id, in_path, out_path):
    # TODO: now only using one molecular model A for the analysis, change it to
    # TODO: use more chains
    in_name = "{}/{}_@_{}.log".format(in_path, pdb_id, subs)
    out_name = "{}/{}-A_{}-A_{}_{}.log".format(out_path, pdb_id, ref_pdb,
                                                 mcsa_id, subs)
    shutil.copy(in_name, out_name)

def _create_substrate(subs, kb_path):
    subs_name = "@substrate,{}".format(subs)
    substrate = {"__init__": [subs, subs], "ref_id": subs_name, "provenance": []}
    properties = {"database": [{"value": "manual", "provenance": []}]}
    properties["cid"] = [{"value": subs, "provenance": [{"source": "manual"}]}]
    properties["query"] = [{"value": subs, "provenance": [{"source": "manual"}]}]
    properties["name"] = [{"value": subs, "provenance": [{"source": "manual"}]}]
    substrate["properties"] = properties
    with open("{}/substrate/{}.json".format(kb_path, subs_name), "w") as out_stream:
        json.dump(substrate, out_stream, indent=2)

def _get_residues_manually(config, jsondb):
    import pandas as pd
    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser(PERMISSIVE=1)
    pdbpath = "{}/{}_c0.pdb".format(config["INPATHS"]["proteins_in_dir"], jsondb["QUERY_PROTEIN"]["PDB"])
    assert os.path.isfile(pdbpath)
    prot = parser.get_structure("prot", pdbpath)
    residues = []
    for mi in jsondb["MAPPING_INFO"]:
        for chain in prot[0].get_chains():
            if not chain.id == mi[0]:
                continue
            for res in prot[0][chain.id].get_residues():
                if res.id[1] == mi[1]:
                    residues.append((mi[1], res.get_resname(), "electrostatic stabiliser"))
    return pd.DataFrame(residues, columns=["resid/chebi id", "PDB code", "ROLE"])

def _get_roles(roles):
    chunks = roles.split("' '")
    roles_list = []
    for chunk in chunks:
        roles_list.append(chunk.strip("[").strip("]").strip("'"))
    return roles_list

def _create_mcsa(mcsa_id, chebi_id, kb_path, config):
    with open(config["INPATHS"]["database_file"], "r") as f:
        jsondb = json.load(f)
    pdb_id = jsondb["QUERY_PROTEIN"]["PDB"]
    
    mcsa_query = "M{:>04d}".format(int(mcsa_id))
    mcsa_entry = {"__init__": [mcsa_id], "ref_id": "@m-csa,{}".format(mcsa_id),
                  "provenance": [{"source": "m-csa"}]}
    pdb_entry = {"value": pdb_id, "provenance": [{"source": "m-csa"}]}
    chebi_entry = {"value": chebi_id, "provenance": [{"source": "m-csa"}]}
    resid_entry = {"provenance": [{"source": "m-csa"}], "value": []}
    
    residues = _get_residues_manually(config, jsondb)
    for index, row in residues.iterrows():
        res = {"position": int(row["resid/chebi id"]),
               "code": row["PDB code"].upper(), "description": "To fill"}
        res["types"] = _get_roles(row["ROLE"])
        resid_entry["value"].append(res)
    properties = {"residues": [resid_entry], "chebi": [chebi_entry],
                  "pdb": [pdb_entry]}
    mcsa_entry["properties"] = properties
    with open("{}/m-csa/{}.json".format(kb_path, "@m-csa,{}".format(mcsa_id)), "w") as out_stream:
        json.dump(mcsa_entry, out_stream, indent=2)

def fill_knowledge_base(kb_path, config, conformers):
    datafile = config["INPATHS"]["database_file"]
    with open(datafile) as in_stream:
        database = json.load(in_stream)
    map_info = database["MAPPING_INFO"]
    active_resids = [int(x[1]) for x in map_info]
    active_resids.sort()
    docked_dir = config["INPATHS"]["docked_dir"]
    subs_in_sdf_dir = config["INPATHS"]["substrates_in_sdf_dir"]
    subs_in_pdbqt_dir = config["INPATHS"]["substrates_in_pdbqt_dir"]
    subs_out_dir = config["OUTPATHS"]["substrates_out_dir"]
    logs_path = config["INPATHS"]["docked_logs"]
    proteins_in_dir = config["INPATHS"]["proteins_in_dir"]
    proteins_out_dir = config["OUTPATHS"]["proteins_out_dir"]
    complexes_path = config["OUTPATHS"]["complexes_dir"]
    for uid, uniprot in database["UNIPROT"].items():
        for pdb_id in uniprot["pdbs"]:
            # Create the mapping json file for each pdb entry in mapping dir
            _create_mapping(pdb_id, database["REF_PDB"], config["EC_INFO"]["mcsa_id"], database, uid, kb_path)
            _create_protein(pdb_id, database["REF_PDB"], config["EC_INFO"]["mcsa_id"], proteins_in_dir, proteins_out_dir)
            
            # Create and fill the complex dir inside knowledge_base, at the same
            # time create complexes dir in data directory and fill it with the
            # docked results.
            for subs in uniprot["substrates"]:
                _create_json_complex(pdb_id, subs, uid, database["REF_PDB"], config["EC_INFO"]["mcsa_id"], kb_path)
                _create_complex(pdb_id, subs, database["REF_PDB"],
                                config["EC_INFO"]["mcsa_id"], docked_dir,
                                complexes_path)
                # _copy_log(pdb_id, subs, database["REF_PDB"],
                #           config["EC_INFO"]["mcsa_id"], logs_path, docked_dir)
        
        # Create and fill the substrate dir
        for subs in uniprot["substrates"]:
            _create_substrate(subs, kb_path)
            sdf_in = "{}/{}.sdf".format(subs_in_sdf_dir, subs)
            sdf_out = "{}/{}.sdf".format(subs_out_dir, subs)
            shutil.copy(sdf_in, sdf_out)
            pdbqt_in = "{}/{}.pdbqt".format(subs_in_pdbqt_dir, subs)
            pdbqt_out = "{}/{}.pdbqt".format(subs_out_dir, subs)
            shutil.copy(pdbqt_in, pdbqt_out)
    
    # Create and fill the m-csa dir
    mcsa_id = config["EC_INFO"]["mcsa_id"]
    mcsa_chebi = config["EC_INFO"]["mcsa_chebi"]
    try:
        config["TASK_INFO"]["manual_dataset"]
    except KeyError:
        raise RuntimeError("Is this not a manual dataset?")
    _create_mcsa(mcsa_id, mcsa_chebi, kb_path, config)
    
    # Create files for the queried file
    for conf in range(conformers):
        query_pdb = "{}_c{}".format(database["QUERY_PROTEIN"]["NAME"], conf)
        _create_mapping(query_pdb, database["REF_PDB"], config["EC_INFO"]["mcsa_id"], database, "None", kb_path, query=True)
        _create_protein(query_pdb, database["REF_PDB"], config["EC_INFO"]["mcsa_id"], proteins_in_dir, proteins_out_dir)

CONFIG_PATTERN = {"EC_INFO": {"mcsa_id":"ID_NUMBER", "mcsa_chebi":ID},
                  "INPATHS": {"fragment_sdf":"../01_docking/cluster_all.sdf",
                              "docked_dir":"docking/outputs",
                              "substrates_in_sdf_dir":"substrates",
                              "substrates_in_pdbqt_dir":"substrates",
                              "docked_logs":"docking/logs",
                              "proteins_in_dir":"receptors",
                              "database_file":None},
                  "OUTPATHS": {"out_path":None,
                               "path_task":None,
                               "knowledge_base":None,
                               "proteins_out_dir":None,
                               "complexes_dir":None,
                               "substrates_out_dir":None},
                  "TASK_INFO": {"chains":"A", "manual_dataset":True}}


def main():
    clusters = 5
    conformers = 1
    for clu in range(clusters):
        clu_path = "cluster_{}".format(clu)
        
        config = deepcopy(CONFIG_PATTERN)
        config["INPATHS"]["database_file"] = os.path.join(clu_path, "subdatabase.json")
        config["OUTPATHS"]["out_path"] = clu_path
        config["OUTPATHS"]["path_task"] = os.path.join(clu_path, "task.json")
        config["OUTPATHS"]["knowledge_base"] = os.path.join(clu_path, "knowledge_base")
        config["OUTPATHS"]["proteins_out_dir"] = os.path.join(clu_path, "proteins")
        config["OUTPATHS"]["complexes_dir"] = os.path.join(clu_path, "complexes")
        config["OUTPATHS"]["substrates_out_dir"] = os.path.join(clu_path, "substrates")
        # Create the file json.task
        make_task_json(config)
        # Create the knowledge_base folder
        knowledge_base_path = config["OUTPATHS"]["knowledge_base"]
        os.makedirs(knowledge_base_path, exist_ok=True)
        os.makedirs("{}/complex".format(knowledge_base_path), exist_ok=True)
        os.makedirs("{}/mapping".format(knowledge_base_path), exist_ok=True)
        os.makedirs("{}/m-csa".format(knowledge_base_path), exist_ok=True)
        os.makedirs("{}/substrate".format(knowledge_base_path), exist_ok=True)
        # Create proteins folder
        os.makedirs(config["OUTPATHS"]["proteins_out_dir"], exist_ok=True)
        # Create the complexes dir
        os.makedirs(config["OUTPATHS"]["complexes_dir"], exist_ok=True)
        os.makedirs(config["OUTPATHS"]["substrates_out_dir"], exist_ok=True)
        
        fill_knowledge_base(knowledge_base_path, config, conformers)

# python3 02_convert_dbs.py
main()
