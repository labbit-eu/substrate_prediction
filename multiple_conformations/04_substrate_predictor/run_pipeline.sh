#!/bin/bash

python3.9 01_build_datasets.py
python3.9 02_substrates_docking_evaluate.py     -c config/CONFIG_FILE.ini > /dev/null
python3.9 03_cluster_sel_residues.py            -c config/CONFIG_FILE.ini > /dev/null
python3.9 04_substrates_cluster_filter.py       -c config/CONFIG_FILE.ini > /dev/null
python3.9 05_substrates_accessibility.py        -c config/CONFIG_FILE.ini > /dev/null
python3.9 06_substrates_filter_create.py        -c config/CONFIG_FILE.ini > /dev/null
python3.9 07_candidates_prepare.py              -c config/CONFIG_FILE.ini > /dev/null
python3.9 08_candidates_docking.py              -c config/CONFIG_FILE.ini > /dev/null
python3.9 09_candidates_docking_evaluate.py     -c config/CONFIG_FILE.ini > /dev/null
python3.9 10_candidates_filter_distances.py     -c config/CONFIG_FILE.ini > /dev/null
python3.9 11_candidates_accessibility.py        -c config/CONFIG_FILE.ini > /dev/null
python3.9 12_candidates_filter_accessibility.py -c config/CONFIG_FILE.ini > /dev/null
python3.9 13_group_filters.py                   -c config/CONFIG_FILE.ini > /dev/null
python3.9 14_check_candidates.py                -c config/CONFIG_FILE.ini
echo "Finished"

