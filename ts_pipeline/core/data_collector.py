#!/usr/bin/env python3

import os, sys, re, argparse, pdb, glob, pprint, json
import ts_pipeline.core.ts_pipeline as tsp
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from openbabel import pybel

def readmol(xyzfile, chrg=0, sanitize=True):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolFromMol2Block(mol.write('mol2'), sanitize=sanitize)
    return mol
    
def remove_trailing_slash(directory):
    if directory != '/' and directory.endswith('/'):
        return directory.rstrip('/')
    return directory

def main():
    parser = argparse.ArgumentParser(
                        prog='TS_pipeline',
                        description='Data scraper for TS_pipeline',)
    parser.add_argument('dirname', nargs='*', help='Directory(-ies) to scrape data from')
    parser.add_argument('-m', '--mode', nargs='?', help='''TS mode (coordinate). Available options:
                        C-C=C-C       - stilbenes;
                        C-C=N-C       - imines;
                        C1C=CCC=C1    - norbornadienes;
                        C1=CCNC=C1    - diarylethenes''', required=True)
    parser.add_argument('-d', '--dump', nargs='?', help='Filename for json dump of final dictionary', 
                        default = f'{os.path.basename(os.getcwd())}')
    parser.add_argument('--individual', action = 'store_true', help='Write individual json dumps for each molecule instead on accumulated one')
    parser.add_argument('--no-pysis', action = 'store_true', help='Collect ORCA (DFT) data only')
    parser.add_argument('--triple-crest', action = 'store_true', help='Separate procedure for triple-crest')
    parser.add_argument('-v', '--verbose', action = 'store_true', help='Verbose printing for debug')
    args = parser.parse_args()
    
    if not args.dirname:
        print('Since no directories provided, trying to scrape from all directories in the current working directory.')
        cwd = os.getcwd()
        directories = [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]
        print(f'The following directories are found and are to be scraped from: {directories}')
        mols = directories
    else:
        mols = args.dirname
        args.individual = True
    if '.json' in args.dump:
        args.dump = os.path.splitext(path)[0]
    mols = list(map(remove_trailing_slash, mols))
    data_collector = {}
    if not args.triple_crest:
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
                if not args.no_pysis:
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
                        
        ####### 3 Checking orca generated data
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
                        case 'C-C=N-C' | 'Greenfield C-C=N-C':
                            dih_value_TS = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *dihedral_nums)
                            dih_value_reac1 = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *dihedral_nums)
                            dih_value_reac2 = tsp.get_dihedral(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *dihedral_nums)
                            ang_value_TS = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_2.xyz', *dihedral_nums[1:])
                            ang_value_reac1 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_4.xyz', *dihedral_nums[1:])
                            ang_value_reac2 = tsp.get_angle(f'{curr_orca_dir}/cregened_conformers_Compound_5.xyz', *dihedral_nums[1:])
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer key value'] = ang_value_TS
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer key value'] = ang_value_reac1
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer key value'] = ang_value_reac2
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca TS conformer dih value'] = dih_value_TS
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac1 conformer dih value'] = dih_value_reac1
                            data_collector[name]['orca conformers properties'][curr_orca_dir]['orca reac2 conformer dih value'] = dih_value_reac2
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
                if args.verbose: print(e)
            except IndexError as e:
                if args.verbose: print(e)
        if args.individual == False:
            print(f'Final length of the data collector {len(data_collector)}, that includes: {mols}')
            with open(f'{args.dump}.json', 'w') as json_file:
                json.dump(data_collector, json_file, indent=4)
        else:
            for key, value in data_collector.items():
                individual_dict = {}
                individual_dict[key] = value
                with open(f'{key}.json', 'w') as json_file:
                    json.dump(individual_dict, json_file, indent=4)
        ####### 4 Processing collected data
    #    conformers_by_orca_count = df['pysis output conformers number'].sum()
    #    conformers_to_orca_count = df['orca input conformers number'].sum()
        conformers_to_orca_total = 0
        conformers_from_orca_total = 0
        failed_orca_count_total = 0
        failed_orca_list_total = []
        energies_dict = {}
        for mol, d in data_collector.items():
            metastable_Es = []
            TS_Es = []
            conformers_to_orca = 0
            conformers_from_orca = 0
            failed_orca_count = 0
            failed_orca_list = []
            for j, conf in enumerate(d['orca output conformers']):
                try:
                    # Checks if properties entry exists for given conf and if all properties dictionary are non-empty
                    if not d['orca conformers properties'][conf] or sum(1 for value in d['orca conformers properties'][conf].values() if value) != len(d['orca conformers properties'][conf].values()):
                        if args.verbose: print(f'Some properties for conformer {conf} are not found, skipping')
                        failed_orca_count += 1
                        failed_orca_list.append(conf)
                        continue
                except KeyError as e:
                    print(f'Conformer {conf} is missing')
                    if args.verbose: print(e)
                    continue
                match args.mode:
                    # Selecting isomer with smaller dihedral as metastable
                    case 'C-C=C-C':
                        if abs(d['orca conformers properties'][conf]['orca reac1 conformer key value']) < abs(d['orca conformers properties'][conf]['orca reac2 conformer key value']):
                            metastable = 'reac1'
                        else:
                            metastable = 'reac2'
                        if args.verbose: print('For conf {conf}, structure with C-C=C-C angle', 
                              d['orca conformers properties'][conf][f'orca {metastable} conformer key value'],
                              'selected as metastable')
                    case 'C-C=N-C' | 'Greenfield C-C=N-C':
                        if abs(d['orca conformers properties'][conf]['orca reac1 conformer dih value']) < abs(d['orca conformers properties'][conf]['orca reac2 conformer dih value']):
                            metastable = 'reac1'
                        else:
                            metastable = 'reac2'
                        if args.verbose: print('For conf {conf}, structure with C=N-C angle', 
                              d['orca conformers properties'][conf][f'orca {metastable} conformer key value'],
                              'and C-C=N-C dihedral',
                              d['orca conformers properties'][conf][f'orca {metastable} conformer dih value'],
                              'selected as metastable')
                    # Selecting isomer with shorter C-C distances as metastable
                    case 'C1C=CCC=C1':
                        if all(d['orca conformers properties'][conf]['orca reac1 conformer key value'][n] < 1.8 for n in [0,1]):
                            metastable = 'reac1'
                        else:
                            metastable = 'reac2'
                        if args.verbose: print('Structure with shorter C-C distances', 
                              d['orca conformers properties'][conf][f'orca {metastable} conformer key value'],
                              'selected as metastable')                    
                    case 'C1=CCNC=C1':
                        if d['orca conformers properties'][conf]['orca reac1 conformer key value'] < d['orca conformers properties'][conf]['orca reac2 conformer key value']:
                            metastable = 'reac1'
                        else:
                            metastable = 'reac2'
                        if args.verbose: print('Structure with shorter C-C distance', 
                              d['orca conformers properties'][conf][f'orca {metastable} conformer key value'],
                              'selected as metastable')                    
                metastable_energy = d['orca conformers properties'][conf][f'orca {metastable} def2-tzvp full energy'] + d['orca conformers properties'][conf][f'orca {metastable} E corr']
                TS_energy =  d['orca conformers properties'][conf]['orca TS def2-tzvp full energy'] + d['orca conformers properties'][conf]['orca TS E corr']
                metastable_Es.append(metastable_energy)
                TS_Es.append(TS_energy)
            failed_orca_list_total.extend(failed_orca_list)
            conformers_to_orca_total += d['orca input conformers number']
            conformers_from_orca_total += len(TS_Es)
            failed_orca_count_total += failed_orca_count
            energies_dict[f"{mol}"] = {}
            energies_dict[f"{mol}"]['Metastable energies'] = metastable_Es
            energies_dict[f"{mol}"]['TS energies'] = TS_Es
        print(f'Conformers generated by ORCA stage: {conformers_to_orca_total}')
        print(f'Successful ORCA conformers: {conformers_from_orca_total}')
        print(f'Probably you should manually check failed ORCA jobs (total {failed_orca_count_total}):')
        print(*[f'{failed_conf}\n' for failed_conf in failed_orca_list_total], end='', sep='')
        if args.verbose: pprint.pprint(energies_dict)
        print(f'Total number of molecules in energies dictionary: {len(energies_dict.keys())}')
        if args.individual == False:
            print(f'Energies are dumped in {args.dump}_energies.json')
            with open(f'{args.dump}_energies.json', 'w') as json_file:
                json.dump(energies_dict, json_file, indent=4)
        else:
            for key, value in energies_dict.items():
                individual_dict = {}
                individual_dict[key] = value
                with open(f'{key}_energies.json', 'w') as json_file:
                    json.dump(individual_dict, json_file, indent=4)
                print(f'Energies for {key} are dumped in {key}_energies.json')
        ######## 5 Williams 
        Gs = {}
        for mol_name, e_dict in energies_dict.items():
            try:
                Gs[mol_name] = tsp.Williams_proc1(mol_name, e_dict, verbose = args.verbose)
            except ValueError as e:
                print(e) 
        if args.verbose: print(*[f'{name: <15} {Geff: <10.2f}' for name, Geff in Gs.items()], sep='\n')
        if args.individual == False:
            with open(f'{args.dump}_Geffs.json', 'w') as json_file:
                json.dump(Gs, json_file, indent=4)
                print(f"Final Williams-treated Geff's are dumped in {args.dump}_Geffs.json")
        else:
            for key, value in Gs.items():
                individual_dict = {}
                individual_dict[key] = value
                with open(f'{key}_Geffs.json', 'w') as json_file:
                    json.dump(individual_dict, json_file, indent=4)
                print(f"Final Williams-treated Geff's for {key} are dumped in {key}_Geffs.json")
    # else is for args.triple_crest == True:
    else:
        energies_dict = {}
        for mol in mols:
            mol_path = f'{mol}/{mol}.xyz'
            match args.mode:
                case 'C-C=C-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    if args.verbose: print(f'TS mode {args.mode} invoked, dihedral angle numbers for {mol}: {dihedral_nums}')
                case 'C-C=N-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    if args.verbose: print(f'TS mode {args.mode} invoked, dihedral angle numbers for {mol}: {dihedral_nums}')
                case 'Greenfield C-C=N-C':
                    dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, 'C-C=N-C', sanitize = False)))
                    if args.verbose: print(f'TS mode {args.mode} invoked, dihedral angle numbers for {mol}: {dihedral_nums}')
                case 'C1C=CCC=C1':
                    cyclohexadiene_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol_path, args.mode, sanitize = False)))
                    if args.verbose: print(f'TS mode {args.mode} invoked, cyclohexadiene_nums numbers for {mol}: {cyclohexadiene_nums}')
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
                    if args.verbose: print(f'TS mode {args.mode} invoked, diarylethene chain numbers numbers for {mol}: {adae_chain_nums}')
                case _:
                    sys.exit('No key TS mode provided (arg -m) or it is not recognized, exiting')
            RSes_paths = glob.glob(f'{mol}/1[01]_*/crest*SP')
            if args.verbose:
                print(f'Overall reactant state conformers for {mol}: {len(RSes_paths)}')
                print(*[f'{path}\n' for path in RSes_paths], sep='')
            Ms = []
            Ss = []
            for conf in RSes_paths:
                curr_conf = glob.glob(f'{conf}/crest_conformers*.xyz')[0]
                match args.mode:
                    case 'C-C=C-C'|'C-C=N-C'|'Greenfield C-C=N-C':
                        key_param = tsp.get_dihedral(curr_conf, *dihedral_nums)
                        with open(f'{conf}/orca_SP.out') as f:
                            orca_text = f.read()
                        E = re.findall(r'(?m)(?s)^FINAL SINGLE POINT ENERGY +(-?\d+\.\d+).*', orca_text)[0]
                        E = float(E)
                        if args.verbose: print(f'Dihedral for {curr_conf}:\n {key_param}, full electronic energy: {E}')
                        if abs(key_param) > 90.0:
                            Ss.append(E)
                        else:
                            Ms.append(E)
                        if args.verbose:
                            print('List of energies for stable state:', Ss)
                            print('List of energies for metastable state:', Ms)
            TSes = []
            TS_paths = glob.glob(f'{mol}/09_*/crest*SP')
            for conf in TS_paths:
                curr_conf = glob.glob(f'{conf}/crest_conformers*.xyz')[0]
                match args.mode:
                    case 'C-C=C-C'|'C-C=N-C'|'Greenfield C-C=N-C':
                        key_param = tsp.get_dihedral(curr_conf, *dihedral_nums)
                        with open(f'{conf}/orca_SP.out') as f:
                            orca_text = f.read()
                        E = re.findall(r'(?m)(?s)^FINAL SINGLE POINT ENERGY +(-?\d+\.\d+).*', orca_text)[0]
                        E = float(E)
                        if args.verbose:
                            print('List of energies for stable state:', Ss)
                            print('List of energies for metastable state:', Ms)
                        TSes.append(E)
                        if args.verbose: print('List of energies for TS:', TSes)
            energies_dict[f'{mol}'] = {'TS energies': TSes,
                                       'Metastable energies': Ms,}
        if args.verbose: 
            print('Final TS and metastable energies:')
            pprint.pprint(energies_dict)
        if args.individual == False:
            print(f'Energies are dumped in {args.dump}_energies.json')
            with open(f'{args.dump}_energies.json', 'w') as json_file:
                json.dump(energies_dict, json_file, indent=4)
        else:
            for key, value in energies_dict.items():
                individual_dict = {}
                individual_dict[key] = value
                with open(f'{key}_energies.json', 'w') as json_file:
                    json.dump(individual_dict, json_file, indent=4)
                print(f'Energies for {key} are dumped in {key}_energies.json')
        # Williams for triple-crest
        Gs = {}
        for mol_name, e_dict in energies_dict.items():
            try:
                Gs[mol_name] = tsp.Williams_proc1(mol_name, e_dict, verbose = args.verbose)
            except ValueError as e:
                print(e) 
        if args.verbose: print(*[f'{name: <15} {Geff: <10.2f}' for name, Geff in Gs.items()], sep='\n')
        if args.individual == False:
            with open(f'{args.dump}_Geffs.json', 'w') as json_file:
                json.dump(Gs, json_file, indent=4)
                print(f"Final Williams-treated Geff's are dumped in {args.dump}_Geffs.json")
        else:
            for key, value in Gs.items():
                individual_dict = {}
                individual_dict[key] = value
                with open(f'{key}_Geffs.json', 'w') as json_file:
                    json.dump(individual_dict, json_file, indent=4)
                print(f"Final Williams-treated Geff's for {key} are dumped in {key}_Geffs.json")

if __name__== '__main__':
    main()