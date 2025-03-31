#!/usr/bin/env python

####### STRUCTURE COMBINING #######
import argparse, os, itertools
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import AllChem, rdMolDescriptors, Draw
f***REMOVED*** rdkit.Chem.rdmolfiles import MolsF***REMOVED***CDXMLFile
f***REMOVED*** collections import defaultdict

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline - fragments combiner',
                    description='Combines molecular fragments f***REMOVED*** two .cdxml documents thorugh U-marked attachment points',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
args = parser.parse_args()

def merge_molecules(mol1, mol2, bond_type=Chem.BondType.SINGLE):
    """
    Merges two molecules by removing the 'U' attachment points and connecting
    the atoms adjacent to them. This version adjusts the neighbor indices based
    on the removal of the U atoms.
    """
    # Convert mol1 to an editable molecule
    combined = Chem.RWMol(mol1)

    # Find attachment points (using 'U' as the placeholder) in both molecules
    u_atoms_1 = [atom for atom in mol1.GetAtoms() if atom.GetSymbol() == "U"]
    u_atoms_2 = [atom for atom in mol2.GetAtoms() if atom.GetSymbol() == "U"]

    if not u_atoms_1 or not u_atoms_2:
        raise ValueError("No U attachment points found in one or both fragments!")

    # For each molecule, record the index of U and its neighbor
    # For mol1:
    u1_idx = u_atoms_1[0].GetIdx()
    neighbor1 = u_atoms_1[0].GetNeighbors()[0]
    n1_idx = neighbor1.GetIdx()
    # Adjust neighbor index: if neighbor comes after U, subtract 1 after removal.
    new_n1_idx = n1_idx - 1 if n1_idx > u1_idx else n1_idx

    # For mol2:
    u2_idx = u_atoms_2[0].GetIdx()
    neighbor2 = u_atoms_2[0].GetNeighbors()[0]
    n2_idx = neighbor2.GetIdx()
    new_n2_idx = n2_idx - 1 if n2_idx > u2_idx else n2_idx

    # Create an editable copy of mol2
    editable_mol2 = Chem.RWMol(mol2)

    # Remove the U atoms f***REMOVED*** both molecules before merging
    editable_mol2.RemoveAtom(u2_idx)  # Remove U f***REMOVED*** mol2
    combined.RemoveAtom(u1_idx)       # Remove U f***REMOVED*** mol1

    # Calculate the offset after removal f***REMOVED*** mol1
    offset = combined.GetNumAtoms()

    # Append atoms f***REMOVED*** mol2 to mol1 (with the adju***REMOVED***d indices f***REMOVED*** mol2)
    for atom in editable_mol2.GetAtoms():
        combined.AddAtom(atom)

    # Append bonds f***REMOVED*** mol2 to mol1 (adjusting indices by the offset)
    for bond in editable_mol2.GetBonds():
        combined.AddBond(
            bond.GetBeginAtomIdx() + offset,
            bond.GetEndAtomIdx() + offset,
            bond.GetBondType()
        )

    # Add a new bond connecting the adju***REMOVED***d neighbor atoms
    combined.AddBond(new_n1_idx, new_n2_idx + offset, bond_type)

    return combined.GetMol()

file_name_count = defaultdict(int)
def get_unique_filename(base_name, ext=".xyz"):
    count = file_name_count[base_name]
    if count > 0:
        unique_name = f"{base_name}_{count}{ext}"
    else:
        unique_name = f"{base_name}{ext}"
    file_name_count[base_name] += 1  # Increment counter
    return unique_name


if len(args.filename) != 2:
    print('List of files submitted:')
    print(f'{file}\n' for file in args.filename)
    print('Please provide exactly TWO fragment files')

molecules_1 = [mol for mol in MolsF***REMOVED***CDXMLFile(args.filename[0]) if mol is not None]
molecules_2 = [mol for mol in MolsF***REMOVED***CDXMLFile(args.filename[1]) if mol is not None]

smiles = []
merged_molecules = []
names = []
#for frag1, frag2 in itertools.product(molecules_1, molecules_2):
for frag1 in molecules_1:
    for frag2 in molecules_2:
        try:
            merged = merge_molecules(frag1, frag2, Chem.BondType.SINGLE)
            merged = Chem.AddHs(merged)
            AllChem.Compute2DCoords(merged)
            merged_molecules.append(merged)
            smiles.append(Chem.MolToSmiles(merged))
        except Exception as e:
            print("Failed to merge a pair of fragments:", e)

for smile in smiles:
    mol = Chem.MolF***REMOVED***Smiles(smile)
    if mol is None:
        print(f"Skipping invalid SMILES: {smile}")
        continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
#    AllChem.MMFFOptimizeMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    file_name = get_unique_filename(formula)
    names.append(file_name)
    xyz_block = Chem.MolToXYZBlock(mol)
    with open(file_name, "w") as f:
        f.write(xyz_block)

    print(f"Saved: {file_name}")

#legends = list(f'{name}\n{smile}' for name, smile in zip(names, smiles))
legends = names
img = Draw.MolsToGridImage(merged_molecules, molsPerRow=3, subImgSize=(300, 200), legends=legends)
img.save('last_generation.png')