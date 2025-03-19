#!/usr/bin/env python3

import sys, os, shutil, glob
import src.TS_pipe as tsp

files = sys.argv[1:]
suspects=[]
for mol in files:
    try:
        key_dihedral = tsp.find_fragment_atoms_old(mol, 'C/C=N\C')
        key_angle = key_dihedral[1:]
        key_angle = list(map(lambda x: x + 1, list(key_angle)))
        key_angle_value = tsp.get_angle(mol, *key_angle)
        key_angle_value = float(key_angle_value)
        if key_angle_value < 160.0:
            print(f'{os.path.dirname(mol): <80} {key_angle}     {key_angle_value: <10.2f} ANGLE IS TOO SMALL FOR INVERSION')
            suspects.append(mol)
        else:
            print(f'{os.path.dirname(mol): <80} {key_angle}     {key_angle_value: <10.2f} angle is ok')
    except ValueError:
        print(f'Apparenty, geometry for {mol} exploded')
        suspects.append(mol)

with open('suspects.txt', 'w') as f:
    for mol in suspects:
        f.write(f'{mol}\n')
#        old_dir = os.path.dirname(mol)
#        new_dir = os.path.basename(os.path.dirname(mol))
#        os.mkdir(new_dir)
#        shutil.copy(mol, new_dir)
#        shutil.copy(f'{old_dir}/TS.yaml', new_dir)
#        
#        to_copy = glob.glob(f'{old_dir}/geom_[12].xyz')
#        print(to_copy)
#        for file in to_copy:
#            shutil.copy(file, new_dir)
