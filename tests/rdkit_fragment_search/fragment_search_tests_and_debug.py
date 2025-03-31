#!/usr/bin/env python3

import os, re, shutil, sys, glob, argparse
f***REMOVED*** openbabel import pybel
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import rdmolfiles, rdDetermineBonds, rdMolTransforms, Draw

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically finds TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
parser.add_argument('-m', '--mode', nargs='?', help='TS mode (coordinate)')
args = parser.parse_args()


def find_fragment_atoms(xyzfile, reference_smiles, chrg=0, sanitize=True):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule f***REMOVED*** an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolF***REMOVED***Mol2Block(mol.write('mol2'), sanitize=sanitize)
    reference_fragment = Chem.MolF***REMOVED***Smiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    if len(matches) > 1:
        raise ValueError("Selected reference fragment (--dihedral) is not unique in the molecule!")
    return matches[0]  # Return the first match    
def find_fragment_atoms_with_hydrogens(xyzfile, reference_smiles, chrg=0, sanitize=False):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule f***REMOVED*** an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolF***REMOVED***Mol2Block(mol.write('mol2'), removeHs = False, sanitize=sanitize)
    reference_fragment = Chem.MolF***REMOVED***Smiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    if len(matches) > 1:
        raise ValueError("Selected reference fragment ({reference_smiles}l) is not unique in the molecule!")
    key_atoms = [mol.GetAtomWithIdx(atom_index) for atom_index in matches[0][1:3]]
    hydrogen_indices = []
    for atom in key_atoms:
        hydrogen_idxs = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"]
        hydrogen_indices.extend(hydrogen_idxs)
    all_atoms = hydrogen_indices + list(matches[0])
    return all_atoms    
    

if __name__ == '__main__':

    sanitize = False
#    reference_smiles = 'N:CCCC:CC'
    reference_smiles = args.mode
    for xyzfile in args.filename:
        print(f'{os.path.basename(xyzfile): <40}')
        mol = next(pybel.readfile('xyz', xyzfile))
        mol = rdmolfiles.MolF***REMOVED***Mol2Block(mol.write('mol2'), sanitize=False)
        count = 0
        smiles_list = ['N:CCCC:CC','N:CC=CC:CC','N=CCCC:CC','N=CC=CC:CC','N:CCCC=CC','N:CC=CC:CC','N=CC=CC=CC', 'N:CC=CC=CC']
        matches=[]
#        breakpoint()
        while not matches:
            reference_fragment = Chem.MolF***REMOVED***Smiles(smiles_list[count])
            matches = mol.GetSubstructMatches(reference_fragment)
            count = count + 1
        print(f'Match found with SMILES {smiles_list[count-1]: <20}: {matches}')
'''        if not matches:
            reference_fragment = Chem.MolF***REMOVED***Smiles('N:CCCC:CC')
            matches = mol.GetSubstructMatches(reference_fragment)
            if not matches
            '''
#        except ValueEror:
#            reference_fragment = Chem.MolF***REMOVED***Smiles('N:CCCC:CC')
#            matches = mol.GetSubstructMatches(reference_fragment)
#        finally:
#            print("Reference fragment not found in the molecule!")
#        if len(matches) > 1:
#            raise ValueError("Selected reference fragment (--dihedral) is not unique in the molecule!")










