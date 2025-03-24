#!/usr/bin/env python3

import glob, os, sys, re
sys.path.append('/proj/scgrp/users/x_***REMOVED***pe/TS_pipeline/src/')
f***REMOVED*** openbabel import pybel
f***REMOVED*** rdkit import Chem
import TS_pipe as tsp
prefix='/proj/scgrp/users/x_***REMOVED***pe/TS_pipeline'

for_fragment_detection = glob.glob(f'{prefix}/data/structures/imines/C*.xyz')
names = list(map(lambda x: re.split(r'\.',os.path.basename(x))[0], for_fragment_detection))
#print(names)
for name in names:
    fragment = tsp.find_fragment_atoms(f'{prefix}/data/structures/imines/{name}.xyz', 'C-C=N-C', sanitize=False)
#    print(fragment)
#    print(name, fragment)
    for conf in glob.glob(f'{prefix}/results/imines/{name}/8_*/cregened_conformers_*/ts_final_geometry.xyz'):
        dih = tsp.get_angle(conf, *[i+1 for i in fragment[1:]])
#        print(dih)
        print(f'{name: <20}{dih}')

