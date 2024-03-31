#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
standardize PDB file by:
1) removing heteroatoms, transforming selenomethionine residues (MSE) to methionines, and keeping the most probable
alternative conformations of atoms  with PyMol script
2) producing biological unit out of the crystal if required PDB header is available - at the moment only 1st biol. unit
is considered

NOTE that a NEW input parameter in the task part of the user-config "knowledge_base.hetero_residues_to_keep" to keep
water molecules or various co-factors in the structure
example: knowledge_base.hetero_residues_to_keep = A427, B382
must be stored as a list of strings knowledge_base.hetero_residues_to_keep = ["A427", "B382"], which is expanded to
((resid 427 and chain A) or (resid 382 and chain B)) in this script
"""

import os
import re

from external import makemultimer
from external import process
from libs import paths
from libs.utils import remove_existing
from libs.configuration import ENV 


def standardize_pdb_file(in_file, out_path, hetero_residues_to_keep=None):
    script_path = os.path.join(paths.temp(), "standardize_pdb_file.py")
    pdb_name = os.path.basename(in_file)
    tmp_file1 = os.path.join(paths.temp(), "tmp_pymol.pdb")
    
    with open(script_path, "w") as output_stream:
        output_stream.write(_prepare_pymol_script(in_file, tmp_file1, hetero_residues_to_keep))

    process.execute([ENV["python_pymol"], script_path], shell=ENV["python_pymol_shell"])

    if not os.path.exists(tmp_file1):
        raise RuntimeError(
            "Failed to standardize PDB file: " + pdb_name + " (no output). "
            "Using the following PyMOl script: " + script_path)

    remove_existing(script_path)
    _copy_header(in_file, tmp_file1, os.path.join(paths.temp(), pdb_name))
    _mm_interface(os.path.join(paths.temp(), pdb_name), out_path)
    remove_existing(tmp_file1)
    remove_existing(os.path.join(paths.temp(), pdb_name))
    

def _prepare_pymol_script(in_file, out_file, hetero_residues_to_keep=None):
    return """import os
import pymol
from pymol import cmd, stored
pymol.pymol_argv = ["pymol", "-qc"]
pymol.finish_launching()
cmd.set("retain_order", 1)

def _expand_residue_to_pymol_selection(residue):
    chain = residue[0]
    resid = residue[1:]
    return "(chain {{:s}} and resid {{:s}})".format(chain, resid)

cmd.load("{:s}","prot", 0,"pdb")
hetero_residues_to_keep = {:s}


# remove solvents and other hetero groups except for those in hetero_residues_to_keep 
if (hetero_residues_to_keep is not None) and len(hetero_residues_to_keep) > 0:
   pymol_selection = "("
   for residue in hetero_residues_to_keep:
       if len(pymol_selection) == 1:
           pymol_selection += _expand_residue_to_pymol_selection(residue)
       else:
           pymol_selection += " or "+_expand_residue_to_pymol_selection(residue)
   pymol_selection += ")"        
else:
    pymol_selection = "none"
cmd.remove("(solvent or organic or inorganic) and not {{:s}}".format(pymol_selection))

#change selenomethionine to methionine
cmd.select("MSEs","resn MSE")
cmd.alter("MSEs", "resn='MET'")
cmd.alter("MSEs", "type='ATOM'")
cmd.alter("MSEs and name SE", "elem='S'")
cmd.alter("MSEs and name SE", "name='SD'")
    
        
#get rid of less probable alternative conformations

#select alternative residues and compute total occupancy & number of atoms in residues
cmd.select("altconf_sele",'not alt ""')
stored.alternative = {{}}
stored.num_atoms = {{}}
cmd.iterate("altconf_sele", "stored.alternative[chain,resi,alt] = 0")
cmd.iterate("altconf_sele", "stored.alternative[chain,resi,alt] +=q")
cmd.iterate("altconf_sele", "stored.num_atoms[chain,resi,alt] = 0")
cmd.iterate("altconf_sele", "stored.num_atoms[chain,resi,alt] += 1")
alt_residues={{}}

#initialize alt_residues
for (chain, resid, alt) in stored.alternative:
    res_sele="chain {{:s}} and resid {{:s}}".format(chain, resid)
    alt_residues[res_sele]={{}}

# assign average occupancy to alt_residues
for (chain, resid, alt) in stored.alternative:
    res_sele="chain {{:s}} and resid {{:s}}".format(chain, resid)
    alt_residues[res_sele][alt]=stored.alternative[chain, resid, alt]/stored.num_atoms[chain, resid, alt]
        
#go through all residues with alternative states
highest={{}}
for res_sele in alt_residues:
    highest["alt"] = "A"
    highest["occupancy"] = 0.0
    #for each residue find the alt_state with highest occupancy 
    for alt in alt_residues[res_sele]:
        if alt_residues[res_sele][alt] > highest["occupancy"]:
            highest["occupancy"] = alt_residues[res_sele][alt]
            highest["alt"] = alt
    #remove all alt_states with lower ocupancy for this reside (or keep A if all states are same)
    cmd.remove(res_sele+" and not (alt ''+{{:s}})".format(highest["alt"]))
    
#reset the pdb info about alternative states
cmd.alter("all", "alt=''" )                
cmd.save("{:s}","prot",-1,"pdb")
cmd.quit()
    """.format(in_file.replace("\\", "/"), str(hetero_residues_to_keep), out_file.replace("\\", "/"))


def _copy_header(original_file, tmp_file, to_file):
    """
    joins PDB header from original PDB file and structural data from standardized PDB file to new PDB file
    """
    header = ""
    with open(tmp_file, "r") as input_stream:
        clean_pdb = input_stream.read()
    
    with open(original_file, "r") as input_stream: 
        for line in input_stream:
            if re.match("^ATOM", line):
                break
            header += line
        full_file_content = header + clean_pdb
    
    with open(to_file, "w") as output_stream:
        output_stream.write(full_file_content)


def _mm_interface(infile, out_path):
    options = {
        "backbone": False,
        "nowater": False,
        "nohetatm": False,
        "renamechains": 1,
        "renumberresidues": 0
    }
    pdblurb = open(infile).read()
    r = makemultimer.PdbReplicator(pdblurb, options)
    outfile_template = os.path.join(out_path, os.path.basename(infile).split('.')[0] + '_mm%s')
    
    for i, bm in enumerate(r.biomolecules):
        outfile = outfile_template % (i+1)
        open(outfile, 'w').write(bm.output(infile))
