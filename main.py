#!/usr/bin/env python3

import os, shutil, argparse
import src.TS_pipe as tsp

### Installation section
pysis_path = '/home/marlen/.local/bin/pysis'
xtb_path = '/usr/local/bin/xtb'
crest_path = '/usr/local/bin/crest3'
orca_path = 'orca.run'
ts_pipe_dir = os.path.dirname(os.path.realpath(__file__))

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically finds TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
args = parser.parse_args()

#%%
if __name__== '__main__':
    mols = args.filename
    mols = list(map(os.path.abspath, mols))

    for mol in mols:
        start_dir = os.getcwd()
        mol_dir = tsp.mkbasedir(mol)
        shutil.copy(mol, mol_dir)
        os.chdir(mol_dir)

####### 1 Optimization
        opt_dir = tsp.mkbasedir(mol, prefix='1_', suffix='_xtb_opt')
        optimized = tsp.xtb_opt(mol, 
                                dirname=opt_dir,  
                                optlev='--optlev extreme')
#                                solvent='--alpb acetonitrile',)
         
######## Detecting dihedral	
        dihedral_nums = tsp.find_fragment_atoms(mol, 'C/C=N\C')
        dihedral_nums = list(map(lambda x: x+1, dihedral_nums))
        opt_dih_value = tsp.get_dihedral(optimized, *dihedral_nums)
        dihedral_line_xtb = ','.join(map(str,dihedral_nums))
        if opt_dih_value < 90:
            dihedral_line_xtb += ',0.0'
            scan_line = '0.0, 180.0, 18'
        else:
            dihedral_line_xtb += ',180.0'
            scan_line = '180.0, 0.0, 18'
        angle_CNC = list(map(str,dihedral_nums[1:]))
        angle_NCC = list(map(str,dihedral_nums[0:3]))
        angle_CNC_line = ','.join(angle_CNC) + ',auto'
        angle_NCC_line = ','.join(angle_NCC) + ',auto'
        
####### 2 Scan dihedral (rotation E to Z or Z to E)
        scan_dih_dir = tsp.mkbasedir(mol, prefix='2_', suffix='_xtb_scan_dih')
        dih_scanned = tsp.xtb_scan_rotation(optimized, 
                                            dirname=scan_dih_dir, 
                                            dihedral=dihedral_line_xtb, 
                                            scan=scan_line,)# solvent='--alpb acetonitrile')
        
####### 3 Second dia***REMOVED***reomer reopt
        m_opt_dir = tsp.mkbasedir(mol, prefix='3_', suffix='_xtb_scan_reopt')
        scan_reoptimized = tsp.xtb_opt(dih_scanned, 
                                       dirname=m_opt_dir,
                                       optlev='--optlev extreme',)
#                                    solvent='--alpb acetonitrile',)

####### 4 Pysis growing string to find TS between E and Z
        for_pysis = tsp.mkbasedir(mol, prefix='4_', suffix='_pysis')
        TS = tsp.pysis_gs(f'{ts_pipe_dir}/templates/TS.yaml',
                          optimized, 
                          scan_reoptimized, 
                          dirname=for_pysis)
        
####### 5 Constrained sampling with CREST
        for_TS_sampling = tsp.mkbasedir(mol, prefix='5_', suffix='_TS_sampling', )
        TS_conformers = tsp.crest_constrained_sampling(TS,
                                                       dirname=for_TS_sampling, 
                                                       angle=[angle_CNC_line, angle_NCC_line],
                                                       optlev='--extreme',)
####################################################### len='--len x3', mdlen='--mdlen x3')
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
                              dirname = for_orca,
                              control_ang = ang_CNC,
                              control_ang_range = [100,140])
        os.chdir(start_dir)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        