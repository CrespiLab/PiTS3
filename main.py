#!/usr/bin/env python3

import os, shutil, argparse
import pdb
import src.TS_pipe as tsp

### Installation section
pysis_path = '/home/marlen/.local/bin/pysis'
xtb_path = '/usr/local/bin/xtb'
crest_path = '/usr/local/bin/crest3'
ts_pipe_dir = os.path.dirname(os.path.realpath(__file__))

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically finds TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
args = parser.parse_args()



if __name__== '__main__':
    mols = args.filename
    mols = list(map(os.path.abspath, mols))

    for mol in mols:
        start_dir = os.getcwd()
        #### SET ignore_existing TO FALSE UNLESS DEBUGGING
        mol_dir = tsp.mkbasedir(mol, ignore_existing=True)
        shutil.copy(mol, mol_dir)
        os.chdir(mol_dir)

# ####### 1 Optimization
#         opt_dir = tsp.mkbasedir(mol, prefix='1_', suffix='_xtb_opt')
#         optimized = tsp.xtb_opt(mol, dirname=opt_dir, solvent='--alpb acetonitrile', 
#                             optlev='--optlev extreme')
         
# ######## Detecting dihedral	
#         print(mol)
#         dihedral_nums = tsp.find_fragment_atoms(mol, 'C/C=N\C')
#         dihedral_nums = list(dihedral_nums)
#         dihedral_nums = list(map(lambda x: x+1, dihedral_nums))
#         dihedral_line = ','.join(map(str,dihedral_nums))
#         dihedral_line_xtb = dihedral_line + ',0.0'
#         angle_CNC_line = ','.join(map(str,dihedral_nums[1:]))
#         angle_CNC_line = angle_CNC_line + ',auto'
        
# ####### 2 Scan dihedral (rotation E to Z or Z to E)
#         scan_dih_dir = tsp.mkbasedir(mol, prefix='2_', suffix='_xtb_scan_dih')
#         dih_scanned = tsp.xtb_scan_rotation(optimized, dirname=scan_dih_dir, 
#                                         dihedral=dihedral_line_xtb, 
#                                         scan='180.0, 0.0, 18', solvent='--alpb acetonitrile')
        
# ####### 3 Second dia***REMOVED***reomer reopt
#         m_opt_dir = tsp.mkbasedir(mol, prefix='3_', suffix='_xtb_scan_reopt')
        scan_reoptimized = tsp.xtb_opt(dih_scanned, dirname=m_opt_dir, solvent='--alpb acetonitrile',
#                                     optlev='--optlev extreme')

#         #4 Pysis growing string to find TS between E and Z
#         for_pysis = tsp.mkbasedir(mol, prefix='4_', suffix='_pysis')
#         TS = tsp.pysis_gs(f'{ts_pipe_dir}/templates/TS.yaml', optimized, scan_reoptimized, dirname=for_pysis)
        
# ####### 5 Constrained sampling with CREST
#         for_TS_sampling = tsp.mkbasedir(mol, prefix='5_', suffix='_TS_sampling', )
#         TS_conformers = tsp.crest_constrained_sampling(TS,
#                                                     dirname=for_TS_sampling, 
#                                                     solvent='--alpb acetonitrile',
#                                                     angle=angle_CNC_line,
#                                                     optlev='--extreme')

        optimized = f'{mol_dir}/1_unsubst_Z_xtb_opt/xtbopt.xyz'
        scan_reoptimized = f'{mol_dir}/3_unsubst_Z_xtb_scan_reopt/xtbopt.xyz'
        TS = f'{mol_dir}/4_unsubst_Z_pysis/ts_final_geometry.xyz'
        TS_conformers = f'{mol_dir}/5_unsubst_Z_TS_sampling/crest_conformers.xyz'

####### 6 Pysis to reoptimize all TS conformers
        for_pysis_TS_reopt = tsp.mkbasedir(mol, prefix='6_', suffix='_pysis_TS_reopt')
        reoptimized_TSes = tsp.pysis_ts_reopt(f'{ts_pipe_dir}/templates/TS_reopt.yaml', TS_conformers, dirname = for_pysis_TS_reopt,concat_breaks=True)
        
####### 7 CREGEN of reoptimized TSes
        for_cregen = tsp.mkbasedir(mol, prefix = '7_', suffix = '_optmized_TS_cregen')
        cregened_ensemble = tsp.cregen(reoptimized_TSes, dirname = for_cregen)
        
####### 8 Pysis IRC
    
####### 8 ORCA wB97x-3c Hess + TSOpt + IRC
        for_orca = tsp.mkbasedir(mol, prefix = '8_', suffix = '_orca')
        
        breakpoint()
                       
        # os.chdir(start_dir)