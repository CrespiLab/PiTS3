#!/usr/bin/env python3

import src.TS_pipe as tsp
import glob, os
f***REMOVED*** openbabel import pybel
f***REMOVED*** rdkit import Chem

xyzfiles = glob.glob('data/structures/stilbenes/C14*.xyz')
for xyzfile in xyzfiles:
    frags = tsp.find_fragment_atoms(xyzfile, 'C-C=C-C')
    dih = tsp.get_dihedral(xyzfile, *[i+1 for i in frags])
    print(f'{os.path.basename(xyzfile): <20} {frags}     {dih}')

