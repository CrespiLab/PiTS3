#!/usr/bin/env python3

import os, sys, re, argparse, pdb, glob, pprint, json
import ts_pipeline.core.ts_pipeline as tsp
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import rdmolfiles
f***REMOVED*** openbabel import pybel

def readmol(xyzfile, chrg=0, sanitize=True):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule f***REMOVED*** an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolF***REMOVED***Mol2Block(mol.write('mol2'), sanitize=sanitize)
    return mol
    
def remove_trailing_slash(directory):
    if directory != '/' and directory.endswith('/'):
        return directory.rstrip('/')
    return directory

def main():
    parser = argparse.ArgumentParser(
                        prog='TS_pipeline',
                        description='Data scraper for TS_pipeline',)
    parser.add_argument('dirname', nargs='?', help='Directory(-ies) to scrape data f***REMOVED***')
    parser.add_argument('-m', '--mode', nargs='?', help='''TS mode (coordinate). Available options:
                        C-C=C-C       - stilbenes
                        C-C=N-C       - imines
                        C1C=CCC=C1    - norbornadienes
                        C1=CCNC=C1    - diarylethenes''', required=True)
    parser.add_argument('-d', '--dump', nargs='?', help='Filename for json dump of final dictionary', 
                        default = f'{os.path.basename(os.getcwd())}_dump.json')
    parser.add_argument('--individual', action = 'store_true', help='Write individual json dumps for each molecule in***REMOVED***ad on accumulated one')
    args = parser.parse_args()
    
    if not args.dirname:
        print('Since no directories provided, trying to scrape f***REMOVED*** all directories in the current working directory.')
        cwd = os.getcwd()
        directories = [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]
        print(f'The following directories are found and are to be scraped f***REMOVED***: {directories}')
        mols = directories
        
    else:
        mols = args.dirname
    mols = list(map(remove_trailing_slash, mols))
    data_collector = {}
    for mol in mols:
        try:
            name=mol
            data_collector[name] = {}
            data_collector[name]['TS mode'] = args.mode
            mol_path = f'{name}/{name}.xyz'
    ####### 0 Detecting key TS node
            match args.mode:
                case 'C-C=C-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    data_collector[name]['key parameter atoms'] = dihedral_nums
                case 'C-C=N-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    data_collector[name]['key parameter atoms'] = dihedral_nums
                case 'Greenfield C-C=N-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, 'C-C=N-C', sanitize = False)))
                    data_collector[name]['key parameter atoms'] = dihedral_nums
                case 'C1C=CCC=C1':
                    cyclohexadiene_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    data_collector[name]['key parameter atoms'] = cyclohexadiene_nums
                case 'C1=CCNC=C1':
                    count = 0
                    smiles_list = ['C1=CCNC=C1', 'C1CCNC=C1']
                    matches=[]
                    while not matches:
                        try:
                            matches = tsp.find_fragment_atoms(mol_path, smiles_list[count], sanitize = False)
                            count += 1
                        except ValueError:
                            print(f'Reference fragment {smiles_list[count-1]} not found, trying {smiles_list[count]}')
                            count += 1
                            continue
                    print(f'Match found with SMILES {smiles_list[count-1]: <20}: {matches}')
                    adae_chain_nums = list(map(lambda x: x+1, matches))
                    data_collector[name]['key parameter atoms'] = adae_chain_nums
                case _:
                    sys.exit('No key TS mode provided (arg -m) or it is not recognized, exiting')
    ####### 1 Checking initial TS geometry        
            first_ts = f'{name}/4_{name}_pysis/ts_final_geometry.xyz'
            match args.mode:
                case 'C-C=C-C':
                    dih_value = tsp.get_dihedral(first_ts, *dihedral_nums)
                    data_collector[name]['TS guess key parameter'] = 'dihedral C-C=C-C'
                    data_collector[name]['TS guess key value'] = dih_value
                case 'C-C=N-C':
                    angle_value = tsp.get_angle(first_ts, *dihedral_nums[1:])
                    data_collector[name]['TS guess key parameter'] = 'angle C=N-C'
                    data_collector[name]['TS guess key value'] = angle_value
                case 'Greenfield C-C=N-C':
                    angle_value = tsp.get_angle(first_ts, *dihedral_nums[1:])
                    data_collector[name]['TS guess key parameter'] = 'angle C=N-C'
                    data_collector[name]['TS guess key value'] = angle_value
                case 'C1C=CCC=C1':
                    distances = [tsp.get_distance(first_ts, *cyclohexadiene_nums[1::4]), tsp.get_distance(first_ts, *cyclohexadiene_nums[2::2])]
                    data_collector[name]['TS guess key parameter'] = 'distances C-C'
                    data_collector[name]['TS guess key value'] = distances
                case 'C1=CCNC=C1':
                    distance = tsp.get_distance(first_ts, *adae_chain_nums[2:4])
                    data_collector[name]['TS guess key parameter'] = 'distances N-C'
                    data_collector[name]['TS guess key value'] = distance

    ####### 2 Checking amount of conformers
            data_collector[name]['pysis input conformers number'] = len(glob.glob(f'{name}/8_*/cregened_conformers*[0-9]'))
            data_collector[name]['pysis output conformers'] = glob.glob(f'{name}/8_*/cregened_conformers*[0-9]/ts_final_geometry.xyz')
            data_collector[name]['pysis output conformers number'] = len(data_collector[name]['pysis output conformers'])
            data_collector[name]['pysis conformers properties'] = {}
            for file in data_collector[name]['pysis output conformers']:
                curr_pysis_conf = file
                data_collector[name]['pysis conformers properties'][curr_pysis_conf] = {}
                match args.mode:
                    case 'C-C=C-C':
                        dih_value = tsp.get_dihedral(curr_pysis_conf, *dihedral_nums)
                        data_collector[name]['pysis conformers properties'][curr_pysis_conf]['TS conformer key value'] = dih_value
                    case 'C-C=N-C':
                        angle_value = tsp.get_angle(curr_pysis_conf, *dihedral_nums[1:])
                        data_collector[name]['pysis conformers properties'][curr_pysis_conf]['TS conformer key value'] = angle_value
                    case 'Greenfield C-C=N-C':
                        angle_value = tsp.get_angle(curr_pysis_conf, *dihedral_nums[1:])
                        data_collector[name]['pysis conformers properties'][curr_pysis_conf]['TS conformer key value'] = angle_value
                    case 'C1C=CCC=C1':
                        distances = [tsp.get_distance(curr_pysis_conf, *cyclohexadiene_nums[1::4]), tsp.get_distance(curr_pysis_conf, *cyclohexadiene_nums[2::2])]
                        data_collector[name]['pysis conformers properties'][curr_pysis_conf]['TS conformer key value'] = distances
                    case 'C1=CCNC=C1':
                        distance = tsp.get_distance(curr_pysis_conf, *adae_chain_nums[2:4])
                        data_collector[name]['pysis conformers properties'][curr_pysis_conf]['TS conformer key value'] = distance
                    
    ####### 3 Checking orca        
            data_collector[name]['orca input conformers number'] = len(glob.glob(f'{name}/9_*/cregened_conformers*reopt'))
            data_collector[name]['orca output conformers'] = glob.glob(f'{name}/9_*/cregened_conformers*reopt')
            data_collector[name]['orca output conformers number'] = len(data_collector[name]['orca output conformers'])
            data_collector[name]['orca conformers properties'] = {}
            for dir in data_collector[name]['orca output conformers']:
                curr_orca_dir = dir
                data_collector[name]['orca conformers properties'][curr_orca_dir] = {}
                with open(f'{curr_orca_dir}/cregened_conformers.out', 'r') as f:
                    text = f.read()
                ts_dG        = re.findall(r'(?m)(?s)COMPOUND JOB  3.*Final Gibbs free energy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  4', text)[0]
                ts_ecorr     = re.findall(r'(?m)(?s)COMPOUND JOB  3.*G-E\(el\) +\.+ +(-?\d+\.\d+).*COMPOUND JOB  4', text)[0]
                ts_def2_e    = re.findall(r'(?m)(?s)COMPOUND JOB  6.*FINAL SINGLE POINT ENERGY +(-?\d+\.\d+).*COMPOUND JOB  7', text)[0]

                reac1_dG     = re.findall(r'(?m)(?s)COMPOUND JOB  4.*Final Gibbs free energy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  5', text)[0]
                reac1_ecorr  = re.findall(r'(?m)(?s)COMPOUND JOB  4.*G-E\(el\) +\.+ +(-?\d+\.\d+).*COMPOUND JOB  5', text)[0]
                reac1_def2_e = re.findall(r'(?m)(?s)COMPOUND JOB  7.*FINAL SINGLE POINT ENERGY +(-?\d+\.\d+).*COMPOUND JOB  8', text)[0]

                reac2_dG     = re.findall(r'(?m)(?s)COMPOUND JOB  5.*Final Gibbs free energy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  6', text)[0]
                reac2_ecorr  = re.findall(r'(?m)(?s)COMPOUND JOB  5.*G-E\(el\) +\.+ +(-?\d+\.\d+).*COMPOUND JOB  6', text)[0]
                reac2_def2_e = re.findall(r'(?m)(?s)COMPOUND JOB  8.*FINAL SINGLE POINT ENERGY +(-?\d+\.\d+)', text)[0]
                
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS Gibbs energy']              =  float(ts_dG)
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS E corr']                    =  float(ts_ecorr)    
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS def2-tzvp full energy']     =  float(ts_def2_e)   

                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 Gibbs energy']           =  float(reac1_dG)    
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 E corr']                 =  float(reac1_ecorr) 
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 def2-tzvp full energy']  =  float(reac1_def2_e)

                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 Gibbs energy']           =  float(reac2_dG)    
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 E corr']                 =  float(reac2_ecorr) 
                data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 def2-tzvp full energy']  =  float(reac2_def2_e)
                if args.mode == 'C1C=CCC=C1':
                    ts_dH        = re.findall(r'(?m)(?s)COMPOUND JOB  3.*Total enthalpy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  4', text)[0]
                    reac1_dH     = re.findall(r'(?m)(?s)COMPOUND JOB  4.*Total enthalpy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  5', text)[0]
                    reac2_dH     = re.findall(r'(?m)(?s)COMPOUND JOB  5.*Total enthalpy +\.+ +(-?\d+\.\d+).*COMPOUND JOB  6', text)[0]
                    data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS enthalpy']              =  float(ts_dH)
                    data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 enthalpy']           =  float(reac1_dH)    
                    data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 enthalpy']           =  float(reac2_dH)    
                match args.mode:
                    case 'C-C=C-C':
                        dih_value_TS = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *dihedral_nums)
                        dih_value_reac1 = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *dihedral_nums)
                        dih_value_reac2 = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *dihedral_nums)
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = dih_value_TS
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = dih_value_reac1
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = dih_value_reac2
                    case 'C-C=N-C':
                        ang_value_TS = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *dihedral_nums[1:])
                        ang_value_reac1 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *dihedral_nums[1:])
                        ang_value_reac2 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *dihedral_nums[1:])
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = ang_value_TS
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = ang_value_reac1
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = ang_value_reac2
                    case 'Greenfield C-C=N-C':
                        ang_value_TS = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *dihedral_nums[1:])
                        ang_value_reac1 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *dihedral_nums[1:])
                        ang_value_reac2 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *dihedral_nums[1:])
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = ang_value_TS
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = ang_value_reac1
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = ang_value_reac2
                    case 'C1C=CCC=C1':
                        dist_value_TS = [tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *cyclohexadiene_nums[1::4]), tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *cyclohexadiene_nums[2::2])]
                        dist_value_reac1 = [tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *cyclohexadiene_nums[1::4]), tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *cyclohexadiene_nums[2::2])]
                        dist_value_reac2 = [tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *cyclohexadiene_nums[1::4]), tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *cyclohexadiene_nums[2::2])]
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = dist_value_TS
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = dist_value_reac1
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = dist_value_reac2
                    case 'C1=CCNC=C1':
                        dist_value_TS = tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *adae_chain_nums[2:4])
                        dist_value_reac1 = tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *adae_chain_nums[2:4])
                        dist_value_reac2 = tsp.get_distance(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *adae_chain_nums[2:4])
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = dist_value_TS
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = dist_value_reac1
                        data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = dist_value_reac2
        except FileNotFoundError as e:
            print(e)
        except IndexError as e:
            print(e)
    if args.individual == False:
        print(f'Final length of the data collector for {mols}: {len(data_collector)}')
        with open(args.dump, 'w') as json_file:
            json.dump(data_collector, json_file, indent=4)
    else:
#        breakpoint()
        for key, value in data_collector.items():
            individual_dict = {}
            individual_dict[key] = value
            with open(f'{key}.json', 'w') as json_file:
                json.dump(individual_dict, json_file, indent=4)
if __name__== '__main__':
    main()