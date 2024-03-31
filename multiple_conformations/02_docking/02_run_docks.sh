#!/bin/bash
# -*- coding: utf-8 -*-

PDBS=(PDB CODES TO USE)

for pdb in ${PDBS[@]}; do
    for dk in `ls ${pdb}/configs`; do
        PATH_TO_AUTODOCK_VINA --config ${pdb}/configs/${dk}
    done
done
