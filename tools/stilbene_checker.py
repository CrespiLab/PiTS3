#!/usr/bin/env python3

import glob, os, sys, re
sys.path.append('/proj/scgrp/users/x_***REMOVED***pe/TS_pipeline/src/')
f***REMOVED*** openbabel import pybel
f***REMOVED*** rdkit import Chem
import TS_pipe as tsp
prefix='/proj/scgrp/users/x_***REMOVED***pe/TS_pipeline'

for_fragment_detection = glob.glob(f'{prefix}/data/structures/stilbenes/C*.xyz')
names = list(map(lambda x: re.split(r'\.',os.path.basename(x))[0], for_fragment_detection))
#print(names)
for name in names:
    fragment = tsp.find_fragment_atoms(f'{prefix}/data/structures/stilbenes/{name}.xyz', 'C-C=C-C', sanitize=False)
#    print(fragment)
#    print(name, fragment)
    for conf in glob.glob(f'{prefix}/results/stilbenes/{name}/4_*/ts_final_geometry.xyz'):
        dih = tsp.get_dihedral(conf, *[i+1 for i in fragment])
#        print(dih)
        print(f'{name: <20}{dih}')

