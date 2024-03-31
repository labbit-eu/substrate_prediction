#!/bin/bash
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#

mkdir -p fpockets
for pdb in `ls pdbs`; do
    echo ${pdb}
    PATH_TO_FPOCKET -m 3.1 -f ./pdbs/${pdb}
    mv ./pdbs/${pdb%.pdb}_out fpockets/${pdb%.pdb}_out
done
