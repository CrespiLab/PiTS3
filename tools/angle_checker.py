#!/usr/bin/env python3

import glob, os, sys, re
sys.path.append('/proj/scgrp/users/x_rompe/TS_pipeline/src/')
from openbabel import pybel
from rdkit import Chem
import TS_pipe as tsp

#for_fragment_detection = glob.glob('data/structures/stilbenes/C*.xyz')
#names = list(map(lambda x: re.split(r'\.',os.path.basename(x))[0], for_fragment_detection))
##print(names)
#for name in names:
#    fragment = tsp.find_fragment_atoms(f'data/structures/stilbenes/{name}.xyz', 'C-C=C-C', sanitize=False)
##    print(name, fragment)
#    for conf in glob.glob(f'results/stilbenes/{name}/8_*/cregened_conformers_*/ts_final_geometry.xyz'):
#        dih = tsp.get_dihedral(conf, *[i+1 for i in fragment])
#        print(f'{name: <20}{dih}')

for_fragment_detection = glob.glob('data/structures/greenfield/C*.xyz')
names = list(map(lambda x: re.split(r'\.',os.path.basename(x))[0], for_fragment_detection))
#print(names)
for name in names:
    fragment = tsp.find_fragment_atoms(f'data/structures/greenfield/{name}.xyz', 'C-C=N-C', sanitize=False)
#    print(name, fragment)
    for conf in glob.glob(f'results/greenfield/problematic/{name}/8_*/cregened_conformers_*/ts_final_geometry.xyz'):
        dih = tsp.get_angle(conf, *[i+1 for i in fragment[1:]])
        print(f'{name: <20}{dih}')

