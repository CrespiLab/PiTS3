#!/usr/bin/env python3

import glob, os, sys, re
sys.path.append('/proj/scgrp/users/x_rompe/TS_pipeline/src/')
from openbabel import pybel
from rdkit import Chem
import TS_pipe as tsp
prefix='/proj/scgrp/users/x_rompe/TS_pipeline'


def check_range_status(values, min_val, max_val):
    in_range = [min_val <= v <= max_val for v in values]
    if all(in_range):
        return True, "all values are in range"
    elif not any(in_range):
        return False, "all values are outside"
    else:
        return False, "some values are outside"


for_fragment_detection = glob.glob(f'{prefix}/data/structures/imines/C*.xyz')
names = list(map(lambda x: re.split(r'\.',os.path.basename(x))[0], for_fragment_detection))
#print(names)
failed = []
flat = []
ok = []
for name in names:
    if len(os.listdir(f'{name}')) < 10:
        print(f'{name: <30} : contains less than 10 directories, skipping')
        failed.append(name)
        continue
    fragment = tsp.find_fragment_atoms(f'{prefix}/data/structures/imines/{name}.xyz', 'C-C=N-C', sanitize=False)
    conformer_angles = []
    for conf in glob.glob(f'{prefix}/results/imines/{name}/8_*/cregened_conformers_*/ts_final_geometry.xyz'):
        dih = tsp.get_angle(conf, *[i+1 for i in fragment[1:]])
        conformer_angles.append(dih)
#        print(f'{name: <20}{dih}')
    print(f'{name: <30} :', check_range_status(conformer_angles, 160.0, 180.0)[1])
#    if not_all_in_range(conformer_angles, 160.0, 180.0):
#        print(f'{name} is fine')
#    else:
#        print(f'{name} generated some wrong TS conformers: ', *[f'{angle}\n' for angle in conformer_angles])

