#!/bin/bash
# -*- coding: utf-8 -*-
#
#  02_run_docks.sh
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


PDBS=(PDB CODES TO USE)

for pdb in ${PDBS[@]}; do
    for dk in `ls ${pdb}/configs`; do
        PATH_TO_AUTODOCK_VINA --config ${pdb}/configs/${dk}
    done
done
