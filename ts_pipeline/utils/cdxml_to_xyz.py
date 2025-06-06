#!/usr/bin/env python

####### STRUCTURE COMBINING #######
import argparse, os, itertools
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Draw
from rdkit.Chem.rdmolfiles import MolsFromCDXMLFile
from collections import defaultdict

file_name_count = defaultdict(int)
def get_unique_filename(base_name, ext=".xyz"):
    count = file_name_count[base_name]
    if count > 0:
        unique_name = f"{base_name}_{count}{ext}"
    else:
        unique_name = f"{base_name}{ext}"
    file_name_count[base_name] += 1  # Increment counter
    return unique_name

def main():
    ### Filename parser
    parser = argparse.ArgumentParser(
                        prog='TS_pipeline - fragments combiner',
                        description='Combines molecular fragments from two .cdxml documents thorugh U-marked attachment points',)

    parser.add_argument('filename', nargs='+', help='XYZ file to process')
    args = parser.parse_args()
    for cdxml_file in args.filename:
        smiles = [Chem.MolToSmiles(mol) for mol in MolsFromCDXMLFile(cdxml_file) if mol is not None]
        names = []
        molecules = []
        for smile in smiles:
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                print(f"Skipping invalid SMILES: {smile}")
                continue
            mol = Chem.AddHs(mol)
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            AllChem.UFFOptimizeMolecule(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            molecules.append(mol)
            file_name = get_unique_filename(formula)
            names.append(file_name)
            molecules
            xyz_block = Chem.MolToXYZBlock(mol)
            with open(file_name, "w") as f:
                f.write(xyz_block)
            print(f"Saved: {file_name}")
            legends = names
            img = Draw.MolsToGridImage(molecules, molsPerRow=3, subImgSize=(300, 200), legends=legends)
            img.save('last_conversion.png')
            
if __name__ == '__main__':
    main()
