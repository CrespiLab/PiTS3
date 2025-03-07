#!/usr/bin/env python3
import os, re, shutil, sys, argparse
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import Draw, rdmolfiles, rdDetermineBonds
import pdb


### Installation section
pysis_path = '/home/marlen/.local/bin/pysis'
xtb_path = '/usr/local/bin/xtb'
crest_path = '/usr/local/bin/crest3'

### Filename parser
parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically find TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
args = parser.parse_args()



#%%
### Functions section
def mkbasedir(path_to_xyzfile, prefix='', suffix='', ignore_existing=False):
    '''
    Creates a directory with the basename of a file
        Arguments:
    filename (str) - name of the file to take basename f***REMOVED***
    suffix (str)   - modify the basename with the suffix (no predefined separator!)
    prefix (str)   - modify the basename with the prefix (no predefined separator!)
        Returns:
    path (str) - path to copied file in the created directory
    '''
    filename = os.path.basename(path_to_xyzfile)
    print(filename)
    split_name = re.split(r'\.', filename)
    dirname = split_name[0]
    if suffix or prefix:
        dirname = prefix + dirname + suffix
    if ignore_existing:
        try:
            os.mkdir(dirname)
        except FileExistsError:
            print(f'Directory {dirname} exists, but still trying to proceed')
    else:
        os.mkdir(dirname)
    path = os.path.abspath(dirname)
    print(f'{dirname} created')
    return path
def safe_dir(file, dirname, rename=''):
    '''
    Moves shell process one layer deeper into the file tree (to prevent overwriting files)
        Arguments:
    file     - file to move
    dir      - dir to move and cd into
        Returns:
    init_path - initial directory
    '''
    init_path = os.getcwd()
    if rename:
        shutil.copy(file, dirname + '/' + rename)
        input_file = os.path.basename(dirname + '/' + rename)
    else:
        shutil.copy(file, dirname)
        input_file = os.path.basename(file)
    os.chdir(dirname)
    return input_file, init_path
def xtb_opt(input_file, dirname='.', model='--gfn2', solvent='', optlev=''):
    '''
    Optimizes molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - directory to run xTB job in
        Returns:
    opt_path - XYZ geometry of optimized molecule
    '''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    os.sy***REMOVED***m(f'{xtb_path} {input_file} --opt {model} {solvent} {optlev} | tee xtb_opt_stdout.log')
    opt_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return opt_path
def xtb_scan_rotation(input_file, dirname='.', model='--gfn2', solvent='',
                      dihedral='0,0,0,0,0', scan='0,0,0', optlev=''):
    '''
    Run a scan along selected dihedral of molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    rotated_path - XYZ geometry of alternative dia***REMOVED***reomer
    '''
    scan_input=f'''$constrain
 force constant=0.05
 dihedral: {dihedral}
$scan
 1: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.sy***REMOVED***m(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} | tee xtb_dih_scan_stdout.log')
    rotated_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return rotated_path
def xtb_scan_inversion(input_file, dirname='.', model='--gfn2', solvent='', angle='0,0,0,0.0', scan='0,0,0', optlev=''):
    '''
    Run a scan along selected angle of molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    inverted_path - XYZ geometry of alternative dia***REMOVED***reomer
    '''
# dihedral: 11,10,8,7,180.0
    scan_input=f'''$constrain
 force constant=5.00
 angle: {angle}
$scan
 1: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.sy***REMOVED***m(f'{xtb_path} {path} --opt --input scan.inp {model} {solvent} {optlev} | tee {path}_xtb_angle_scan.log')
    inverted_path = dir + '/xtbopt.xyz'
    os.chdir(initial_path)
    return inverted_path
def pwd(prnt=False):
    '''
    Just python-wrapped pwd
    '''
    os.sy***REMOVED***m('pwd')
def pysis_gs(input_yaml, xyz_1, xyz_2, dirname='.'):
    '''
    Calls pysis (pysisyphus) growing string method with .yaml input and geometries xyz_1 and xyz_2
        Arguments:
    input_yaml   - YAML specification of the job (details at https://pysisyphus.readthedocs.io/en/latest/index.html)
    xyz_1        - first geometry (PES minimum)
    xyz_2        - second geometry (PES minimum)
    '''
    if not input_yaml: print("Please provide input .yaml specification for pysis")
    if not (xyz_1 and xyz_2): 
        print("You need exactly two geometries to continue, exiting")
        sys.exit()
    input_yaml, initial_path = safe_dir(input_yaml, dirname)    
    xyz_1 = safe_dir(xyz_1, dirname, rename='geom_1.xyz')[0]
    xyz_2 = safe_dir(xyz_2, dirname, rename='geom_2.xyz')[0]
    os.chdir(dirname)
    os.sy***REMOVED***m(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
    os.chdir(initial_path)
    TS_path = f'{dirname}/ts_final_geometry.xyz'
    return TS_path
def pysis_ts_reopt(input_yaml, xyz, dirname='.'):
    '''
    Calls pysis (pysisyphus) hessian-based method with .yaml input and geometry xyz1
        Arguments:
    input_yaml   - YAML specification of the job (details at https://pysisyphus.readthedocs.io/en/latest/index.html)
    xyz          - geometry of a TS or multi-xyz
    '''
    if not input_yaml: print("Please provide input .yaml specification for pysis")
    if not xyz: 
        print("You need a TS guess(es) to continue, exiting")
        sys.exit()
    input_yaml, init_path_top = safe_dir(input_yaml, dirname)    
    shutil.copy(xyz, dirname)    
    os.chdir(dirname)
    conformers_file = os.path.basename(xyz)
    TS_conformers = split_geoms(conformers_file)[0]
    reoptimized_TSes = []
    for conformer in TS_conformers:
        curr_conf_dir = mkbasedir(conformer)
        curr_TS_conformer, init_path = safe_dir(conformer, curr_conf_dir)
        os.chdir(curr_conf_dir)
        shutil.copy(f'../{input_yaml}', '.')
        replace_in_file(f'{input_yaml}', 'xyzfile.xyz', f'{curr_TS_conformer}')
        os.sy***REMOVED***m(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
        reoptimized_TSes.append(curr_TS_conformer)
        os.chdir(init_path)
    os.chdir(init_path_top)
    # TS_path = f'{dirname}/ts_final_geometry.xyz'
    return reoptimized_TSes
def find_fragment_atoms(xyzfile, reference_smiles):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule f***REMOVED*** an XYZ file.

    Fixes:
    - Uses RDKit's DetermineConnectivity() to assign bonds to XYZ molecules.
    - Sanitizes the molecule to ensure valence is properly assigned.
    
    :param mol: RDKit molecule
    :param reference_fragment: RDKit molecule representing the common fragment
    :return: List of atom indices that match the fragment
    """
    mol = rdmolfiles.MolF***REMOVED***XYZFile(xyzfile)
    rdDetermineBonds.DetermineConnectivity(mol)
    Chem.rdDetermineBonds.DetermineBondOrders(mol)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    reference_fragment = Chem.MolF***REMOVED***Smiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    return matches[0]  # Return the first match    
def crest_constrained_sampling(input_file,
                               dirname='.',
                               model='--gfn2', 
                               solvent='', 
                               angle='0,0,0,0.0', 
                               optlev='', 
                               cinp='constraints.inp'):
    '''
    Runs CREST conformational sampling with constrain (predefined in the function).
    '''
    input_file, initial_path = safe_dir(input_file, dirname)    
    os.chdir(dirname)
    constraint_input=f'''$constrain
 force constant=5.00
 angle: {angle}
 dihedral: 5,7,8,10
$end'''
    with open('constraints.inp','w') as file:
        file.write(constraint_input)
    crest_line=f'{crest_path} {input_file} --cinp {cinp} {optlev} {model} {solvent} | tee xtb_TS_conf_sampling_stdout.log'
    os.sy***REMOVED***m(crest_line)
    TS_conformers_path = f'{dirname}/crest_conformers.xyz'
    os.chdir(initial_path)
    return TS_conformers_path
geometry_regex = re.compile(r'(?m)\s*\d+\s*\s+.+\s+(?:\w+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*)+')
def split_geoms(xyzfile, regex = geometry_regex, write_geometries = True):
    '''
    Splits multixyz file into set of single geometries and writes into separate files
        Arguments:
    xyzfile          - multixyz file (usually .xyz or .trj)
    regex            - regex to detect geometry inside xyz file
        Returns:
    geometries       - list of geometry strings
    filenames        - corresponding split files
    '''
    basename = os.path.basename(xyzfile)
    basename = re.split(r'\.', basename)[0]
    with open(xyzfile, 'r') as file:
        data = file.read()
    geometries = re.findall(regex, data)
    filenames = []
    if write_geometries:
        for i, geom in enumerate(geometries):
            num_digits = len(str(len(geometries)))
            filename = f'{basename}_{i:0{num_digits}}.xyz'
            filenames.append(filename)
            with open(filename, 'w') as file:
                file.write(geom)
    return filenames, geometries
def replace_in_file(file, find, replace):
    with open(file) as r:
        text = r.read().replace(find, replace)
    with open(file, "w") as w:
        w.write(text)

#%%

#######################################################################
#######################################################################
#######################################################################

if __name__== '__main__':
    mols=args.filename

    for mol in mols:
        start_dir = os.getcwd()
        #### SET ignore_existing TO FALSE UNLESS DEBUGGING
        mol_dir = mkbasedir(mol, ignore_existing=True)
        shutil.copy(mol, mol_dir)
        os.chdir(mol_dir)

        # #1 Optimization
        # opt_dir = mkbasedir(mol, prefix='1_', suffix='_xtb_opt')
        # optimized = xtb_opt(mol, dirname=opt_dir, solvent='--alpb acetonitrile', 
        #                     optlev='--optlev extreme')
         
        # Detecting dihedral	
        print(mol)
        dihedral_nums = find_fragment_atoms(mol, 'FC=NC')
        dihedral_nums = list(dihedral_nums)
        dihedral_nums = list(map(lambda x: x+1, dihedral_nums))
        dihedral_line = ','.join(map(str,dihedral_nums))
        dihedral_line_xtb = dihedral_line + ',0.0'
        angle_CNC_line = ','.join(map(str,dihedral_nums[1:]))
        angle_CNC_line = angle_CNC_line + ',auto'
        
        # #2 Scan dihedral (rotation E to Z or Z to E)
        # scan_dih_dir = mkbasedir(mol, prefix='2_', suffix='_xtb_scan_dih')
        # dih_scanned = xtb_scan_rotation(optimized, dirname=scan_dih_dir, 
        #                                 dihedral=dihedral_line_xtb, 
        #                                 scan='0.0, 180.0, 18', solvent='--alpb acetonitrile')
        
        # #3 Second dia***REMOVED***reomer reopt
        # m_opt_dir = mkbasedir(mol, prefix='3_', suffix='_xtb_scan_reopt')
        # scan_reoptimized = xtb_opt(dih_scanned, dirname=m_opt_dir, solvent='--alpb acetonitrile',
        #                             optlev='--optlev extreme')

        optimized = f'{mol_dir}/1_unsubst_Z_xtb_opt/xtbopt.xyz'
        scan_reoptimized = f'{mol_dir}/3_unsubst_Z_xtb_scan_reopt/xtbopt.xyz'
        TS = f'{mol_dir}/4_unsubst_Z_pysis/ts_final_geometry.xyz'
        TS_conformers = f'{mol_dir}/5_unsubst_Z_TS_sampling/crest_conformers.xyz'

        # #4 Pysis growing string to find TS between E and Z
        # for_pysis = mkbasedir(mol, prefix='4_', suffix='_pysis')
        # TS = pysis_gs('../TS.yaml', optimized, scan_reoptimized, dirname=for_pysis)
        
        # #5 Constrained sampling with CREST
        # for_TS_sampling = mkbasedir(mol, prefix='5_', suffix='_TS_sampling', )
        # TS_conformers = crest_constrained_sampling(TS,
        #                                             dirname=for_TS_sampling, 
        #                                             solvent='--alpb acetonitrile',
        #                                             angle=angle_CNC_line,
        #                                             optlev='--extreme')
        
        #6 Pysis to reoptimize all TS conformers
        for_pysis_TS_reopt = mkbasedir(mol, prefix='6_', suffix='_pysis_TS_reopt')
        reoptimized_TSes = pysis_ts_reopt('/mnt/c/Users/***REMOVED***pe672/Documents/GitHub/TS_pipeline/TS_reopt.yaml', TS_conformers, dirname = for_pysis_TS_reopt)
        
        # TS = pysis_ts_reopt('../TS_reopt.yaml', , scan_reoptimized, dirname=for_pysis)
                
        # os.chdir(start_dir)