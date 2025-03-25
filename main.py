#!/usr/bin/env python3

import os, shutil, argparse, sys
import src.TS_pipe as tsp
import pdb

### Installation section
pysis_path = 'pysis'
xtb_path = 'xtb'
crest_path = 'crest3'
orca_path = 'orca.run'
ts_pipe_dir = os.path.dirname(os.path.realpath(__file__))

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically finds TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
parser.add_argument('-m', '--mode', nargs='?', help='TS mode (coordinate)')
args = parser.parse_args()

#%%
if __name__== '__main__':
    mols = args.filename
    mols = list(map(os.path.abspath, mols))


    for mol in mols:

####### 0 Detecting key TS node
        match args.mode:
            case 'C-C=C-C' | 'C-C=N-C':
                dihedral_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol, args.mode, sanitize = False)))
                opt_dih_value = tsp.get_dihedral(optimized, *dihedral_nums)
                dihedral_line_xtb = ','.join(map(str,dihedral_nums))
                if abs(opt_dih_value) < 90:
                    dihedral_line_xtb += ',0.0'
                    scan_line = '0.0, 180.0, 18'
                else:
                    dihedral_line_xtb += ',180.0'
                    scan_line = '180.0, 0.0, 18'
                angle_CNC_nums = dihedral_nums[1:]
                constraint = tsp.find_fragment_atoms_with_hydrogens(mol, args.mode, sanitize=False)
            case 'C1C=CCC=C1':
                cyclohexadiene_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol, args.mode, sanitize = False)))
                distance1, distance2 = cyclohexadiene_nums[1::4], cyclohexadiene_nums[2::2]
                print(f'''Reference smiles {args.mode} found, atom numbers: {cyclohexadiene_nums}, concerted scan will be done between atom pairs {distance1} and {distance2}''')
                constraint = tsp.find_fragment_atoms_with_hydrogens(mol, args.mode, sanitize=False)
            case 'N=CC=CC=C':
                count = 0
                smiles_list = ['N:CC=CC=CC', 'N:CC=CC:CC', 'N:CCCC:CC', 'N=CCCC:CC','N=CC=CC:CC','N:CCCC=CC','N:CC=CC:CC','N=CC=CC=CC',]
                matches=[]
                while not matches:
                    try:
                        matches = tsp.find_fragment_atoms(mol, smiles_list[count], sanitize = False)
                        count += 1
                    except ValueError:
                        print(f'Reference fragment {smiles_list[count-1]} not found, trying {smiles_list[count]}')
                        count += 1
                        continue
                print(f'Match found with SMILES {smiles_list[count-1]: <20}: {matches}')
                adae_chain_nums = list(map(lambda x: x+1, matches))
                distance1 = [adae_chain_nums[0], adae_chain_nums[-2]]
                starting_distance = tsp.get_distance(mol, adae_chain_nums[0], adae_chain_nums[-2])
                print(f'''Reference smiles {args.mode} found, atom numbers: {adae_chain_nums}, concerted scan will be done between atom pairs {distance1}''')
                constraint = tsp.find_fragment_atoms_with_hydrogens(mol, smiles_list[count-1], sanitize=False)
            case _:
                try:
                    fragment_nums = list(map(lambda x: x+1, tsp.find_fragment_atoms(mol, args.mode, sanitize = False)))
                    print(fragment_nums)
                except ValueError:
                    print(f'Reque***REMOVED***d SMILES {args.mode} not found in the molecule, exiting')
                    sys.exit('No key TS mode provided (arg -m) or it is not recognized, exiting')

####### Saving starting dir
        start_dir = os.getcwd()
        mol_dir = tsp.mkbasedir(mol, ignore_existing=True)
        shutil.copy(mol, mol_dir)
        os.chdir(mol_dir)
####### 1 Optimization
        opt_dir = tsp.mkbasedir(mol, prefix='1_', suffix='_xtb_opt')
        optimized = tsp.xtb_opt(mol, 
                                dirname=opt_dir,  
                                optlev='--optlev extreme',
                                etemp='--etemp 1500')
##                                solvent='--alpb acetonitrile',)
         

####### 2 Scan
        scan_dir = tsp.mkbasedir(mol, prefix='2_', suffix='_xtb_scan')
        match args.mode:
            case 'C-C=C-C' | 'C-C=N-C':
                scanned = tsp.xtb_scan_rotation(optimized, 
                                                dirname=scan_dir, 
                                                dihedral=dihedral_line_xtb, 
                                                scan=scan_line,
                                                etemp='--etemp 1500',)# solvent='--alpb acetonitrile')
                        
            case 'C1C=CCC=C1':
                scanned = tsp.xtb_scan_nbd(optimized,
                                          dirname=scan_dir,
                                          distance1 = distance1,
                                          distance2 = distance2,
                                          scan = '2.5, 1.5, 50',
                                          etemp='--etemp 1500',)
            case 'N=CC=CC=C':
                scanned = tsp.xtb_scan_adae(optimized,
                                           dirname=scan_dir,
                                           distance1 = distance1,
                                           scan = f'{starting_distance}, 1.5, 50',
                                           etemp='--etemp 1500',)
        
####### 3 Second dia***REMOVED***reomer reopt
        m_opt_dir = tsp.mkbasedir(mol, prefix='3_', suffix='_xtb_scan_reopt')
        scan_reoptimized = tsp.xtb_opt(scanned, 
                                       dirname=m_opt_dir,
                                       optlev='--optlev extreme',
                                       etemp='--etemp 1500',)
#                                    solvent='--alpb acetonitrile',)

####### 4 Pysis growing string to find TS between E and Z
        for_pysis = tsp.mkbasedir(mol, prefix='4_', suffix='_pysis')
        TS = tsp.pysis_gs(f'{ts_pipe_dir}/templates/TS_10_nodes.yaml',
                          optimized, 
                          scan_reoptimized, 
                          dirname=for_pysis)
'''
#        TS = '/proj/scgrp/users/x_***REMOVED***pe/TS_pipeline/results/new_imines/C13H11N/4_C13H11N_pysis/ts_final_geometry.xyz'
####### 5 Constrained sampling with CREST
        for_TS_sampling = tsp.mkbasedir(mol, prefix='5_', suffix='_TS_sampling', )
        TS_conformers = tsp.crest_constrained_sampling(TS,
                                                       dirname=for_TS_sampling, 
                                                       constrain_atoms=constraint,
                                                       optlev='--extreme',)
######################################################  dlen='--len x3', mdlen='--mdlen x3')
                                                      # solvent='--alpb acetonitrile',
####### 6 Pysis to reoptimize all TS conformers
        for_pysis_TS_reopt = tsp.mkbasedir(mol, prefix='6_', suffix='_pysis_TS_reopt')
        reoptimized_TSes = tsp.pysis_ts_reopt(f'{ts_pipe_dir}/templates/TS_reopt.yaml', 
                                              TS_conformers, 
                                              dirname = for_pysis_TS_reopt,
                                              concat_breaks=True)
        
####### 7 CREGEN of reoptimized TSes
        for_cregen = tsp.mkbasedir(mol, prefix = '7_', suffix = '_optmized_TS_cregen')
        cregened_ensemble = tsp.cregen(reoptimized_TSes, dirname = for_cregen)
        
####### 8 Pysis IRC
        for_pysis_irc = tsp.mkbasedir(mol, prefix = '8_', suffix = '_filtered_TS_irc')
        cregened_ensemble = shutil.copy(cregened_ensemble, f'{for_pysis_irc}/cregened_conformers.xyz')
        ircs = tsp.pysis_ts_irc(f'{ts_pipe_dir}/templates/TS_reopt_irc.yaml',
                                  xyz = cregened_ensemble,
                                  dirname = for_pysis_irc,)

####### 9 ORCA wB97x-3c Hess + TSOpt + IRC
        for_orca = tsp.mkbasedir(mol, prefix = '9_', suffix = '_orca')
        tsp.orca_three_points(ircs,
                              orca_template = f'{ts_pipe_dir}/templates/orca_three_points.inp',
                              dirname = for_orca,)
#                              control_ang = angle_CNC_nums,
#                              control_ang_range = [100,140])
        os.chdir(start_dir)
        
'''     