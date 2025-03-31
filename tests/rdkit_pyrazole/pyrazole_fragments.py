#!/usr/bin/env python3

import os, re, shutil, sys, pdb, glob
f***REMOVED*** openbabel import pybel
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import rdmolfiles, rdDetermineBonds, rdMolTransforms

mols = glob.glob('greenfield/*.xyz')


xyzfile = mols[0]
reference_smiles = 'C-C=N-C'
mol = next(pybel.readfile('xyz', xyzfile))
mol = rdmolfiles.MolF***REMOVED***Mol2Block(mol.write('mol2'))
breakpoint()
reference_fragment = Chem.MolF***REMOVED***Smiles(reference_smiles)
matches = mol.GetSubstructMatches(reference_fragment)
if not matches:
    raise ValueError("Reference fragment not found in the molecule!")
if len(matches) > 1:
    raise ValueError("Selected reference fragment (--dihedral) is not unique in the molecule!")

    




