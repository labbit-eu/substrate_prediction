#!/bin/bash
#
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#

RESIDS=(45 132 166)
export LRIP_SAMPLING="PATH_TO_LRIP_SAMPLING"

mkdir -p pdbs
for res in ${RESIDS[@]}; do
    cat<<EOF>make_pdbs.cppin
parm ${LRIP_SAMPLING}/sim.top [pref]
trajin ${LRIP_SAMPLING}/298k-rip-310k/${res}/md.trj
reference ${LRIP_SAMPLING}/sim.pdb [simref] parm [pref]
rmsd :1-328@N,CA,C,O,H reference [simref]
trajout pdbs/A0LKP2_${res}.pdb pdb multi
go
quit
EOF
    cpptraj -i make_pdbs.cppin
    cat<<EOF>change_names.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright 2023 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>

import os
pdb_files = os.listdir("./pdbs")
for f in pdb_files:
    if ".pdb." in f:
        fname = f.split(".", 3)
        new_name = fname[0] + ".{:0>3d}.pdb".format(int(fname[2]))
        oldf = os.path.join(".", "pdbs", f)
        newf = os.path.join(".", "pdbs", new_name)
        os.rename(oldf, newf)
EOF
    python3 change_names.py
done
rm make_pdbs.cppin change_names.py
