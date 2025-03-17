#!/usr/bin/env python3

import os, re, shutil, sys
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import rdmolfiles, rdDetermineBonds



### Installation section
pysis_path = '/home/marlen/.local/bin/pysis'
xtb_path = '/usr/local/bin/xtb'
crest_path = '/usr/local/bin/crest3'
orca_path = 'orca.run'

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
def ***REMOVED***p():
    def number_generator():
        num = 1
        while True:
            yield num
            num += 1
    gen = number_generator()  # Create the generator instance
    def next_number():
        print('Step ', next(gen), ' of TS_pipe script')  # Automatically print the next number
    return next_number  # Return the function that prints the next number
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
def xtb_scan_inversion(input_file, dirname='.', model='--gfn2', solvent='',
                       angle='0,0,0,0.0', scan='0,0,0', optlev=''):
    '''
    INVERSION FOR MORE THAN 180 DEGREES IS PROBLEMATIC
    Run a scan along selected angle of molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    inverted_path - XYZ geometry of alternative dia***REMOVED***reomer
    '''
    scan_input=f'''$constrain
 force constant=5.00
 angle: {angle}
$scan
 1: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.sy***REMOVED***m(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} | tee xtb_angle_scan.log')
    inverted_path = os.path.abspath('xtbopt.xyz')
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
def pysis_ts_reopt(input_yaml, xyz, dirname = '.', concat_breaks=False, irc=False):
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
    if irc == False:
        reoptimized_TSes = []
    else:
        reoptimized_TSes = {}
    for conformer in TS_conformers:
        curr_conf_dir = mkbasedir(conformer)
        curr_TS_conformer, init_path = safe_dir(conformer, curr_conf_dir)
        os.chdir(curr_conf_dir)
        shutil.copy(f'../{input_yaml}', '.')
        replace_in_file(f'{input_yaml}', 'xyzfile.xyz', f'{curr_TS_conformer}')
        if irc == True:
            replace_in_file(f'{input_yaml}', 'forward: False', 'forward: True')
            replace_in_file(f'{input_yaml}', 'backward: False', '''backward: True\nendopt:''')
        os.sy***REMOVED***m(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
        reoptimized_TS_conformer = re.split('\.', curr_TS_conformer)[0] + '_reopt.xyz'
        # breakpoint()
        if irc == False:
            reoptimized_TS_conformer = shutil.copy('ts_final_geometry.xyz', f'{reoptimized_TS_conformer}')
            reoptimized_TSes.append(os.path.abspath(reoptimized_TS_conformer))
        else:
            reoptimized_TSes[f'{reoptimized_TS_conformer}'] = {'backward' : os.path.abspath('backward_end_opt.xyz'),
                                                               'TS' : os.path.abspath('ts_final_geometry.xyz'),
                                                               'forward' : os.path.abspath('forward_end_opt.xyz')}
        os.chdir(init_path)
    if irc == False:
        reoptimized_TSes = concatenate(f'{os.path.basename(os.getcwd())}_reoptimized_TSes.xyz', reoptimized_TSes, add_line_break=concat_breaks)
        reoptimized_TSes = os.path.abspath(reoptimized_TSes)
    os.chdir(init_path_top)
    return reoptimized_TSes
def pysis_ts_irc(input_yaml, xyz, dirname = '.'):
    '''
    Calls pysis (pysisyphus) hessian-based method with .yaml input and geometry xyz, then runs IRC with endopt
        Arguments:
    input_yaml   - YAML specification of the job (details at https://pysisyphus.readthedocs.io/en/latest/index.html)
    xyz          - geometry of a TS or multi-xyz
    '''
    if not input_yaml: print("Please provide input .yaml specification for pysis")
    if not xyz: 
        print("You need a TS guess(es) to continue, exiting")
        sys.exit()
    input_yaml, init_path_top = safe_dir(input_yaml, dirname)    
    try:
        shutil.copy(xyz, dirname)
    except shutil.SameFileError:
        print(f'File {xyz} already exists and is the same, but continue')
    os.chdir(dirname)
    conformers_file = os.path.basename(xyz)
    TS_conformers = split_geoms(conformers_file)[0]
    reoptimized_TSes = {}
    reoptimized_TSes['basename'] = re.split(r'\.', os.path.basename(conformers_file))[0]
    reoptimized_TSes['conformers'] = {}
    for conformer in TS_conformers:
        curr_conf_dir = mkbasedir(conformer)
        curr_TS_conformer, init_path = safe_dir(conformer, curr_conf_dir)
        os.chdir(curr_conf_dir)
        shutil.copy(f'../{input_yaml}', '.')
        replace_in_file(f'{input_yaml}', 'xyzfile.xyz', f'{curr_TS_conformer}')
        os.sy***REMOVED***m(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
        reoptimized_TS_conformer = re.split('\.', curr_TS_conformer)[0] + '_reopt.xyz'
        # breakpoint()
        reoptimized_TSes['conformers'][f'{reoptimized_TS_conformer}'] = {'backward' : os.path.abspath('backward_end_opt.xyz'),
                                                               'TS' : os.path.abspath('ts_final_geometry.xyz'),
                                                               'forward' : os.path.abspath('forward_end_opt.xyz')}
        os.chdir(init_path)
    os.chdir(init_path_top)
    return reoptimized_TSes
def concatenate(filename, list_of_files, add_line_break=False):
    '''
    Concatenates text files (optionally adds an empty line in the end of each chunk)
    
    Parameters
    ----------
    filename : string
        Target file for concatenation
    list_of_files : list of string
        List of files for concatenation

    Returns
    -------
    Path to the output file with contatenated text
    '''
    with open(filename, 'w') as file:
        for chunk in list_of_files:
            with open(chunk, 'r') as chunk:
                to_append = chunk.read()
                file.write(to_append)
                if add_line_break: file.write('\n')
    return filename
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
def cregen(ensemble_file, extra_params='', dirname = '.', sorted_ensemble=False):
    '''
    Runs standalone CREST CREGEN procedure on an ensemble file

    Parameters
    ----------
    ensemble_file : str
        Ensemble file
    extra_params : str
        Additional parameters for CREGEN run (e.g. --ewin, --ethr; see CREST documentation https://crest-lab.github.io/crest-docs/page/documentation/keywords.html#ensemble-sorting-options)
    sorted : bool
        If True, returns .xyz.sorted file, otherwise just filtered crest_ensemble.xyz 
    Returns
    -------
    Path to CREGENed file.

    '''
    ensemble_file, init_path = safe_dir(ensemble_file, dirname)
    os.sy***REMOVED***m(f'{crest_path} {ensemble_file} --cregen {ensemble_file} {extra_params} | tee cregen_stdout.log')
    basename = os.path.basename(ensemble_file)
    basename = re.split(r'\.', basename)[0]
    if sorted_ensemble == True:
        resulting_ensemble = os.path.abspath(f'{basename}.xyz.sorted')
    else:
        resulting_ensemble = os.path.abspath('crest_ensemble.xyz')
    os.chdir(init_path)
    return resulting_ensemble
def orca_three_points(irc_dict, orca_template = 'orca_three_points.inp', dirname = '.'):
    '''
    Runs ORCA compound job to reoptimize TS and IRC endpoints. DOES NOT DO DFT IRC!

    Parameters
    ----------
    irc_dict : str
        Dictionary f***REMOVED*** pysis_ts_irc, containing paths to TS and both endpoins for each TS conformer.
    orca_template : str
        ORCA compound job template. Do not change if you do not know how to set compound job!
    dirname : TYPE, optional
        Directory in which to run the job. The default is '.'.

    Returns
    -------
    Path to orca output.

    '''
    init_path_top = os.getcwd()
    os.chdir(dirname)
    basename = irc_dict['basename']
    orca_template = shutil.copy(orca_template, f'{basename}.inp')
    for conformer, xyzs in irc_dict['conformers'].items():
        curr_conf_dir = re.split(r'\.', conformer)[0]
        init_path = os.getcwd()
        os.mkdir(curr_conf_dir)
        os.chdir(curr_conf_dir)
        for value in xyzs.values():
            shutil.copy(value, '.')
        curr_orca_template = shutil.copy(f'../{orca_template}', '.')
        replace_in_file(curr_orca_template, '%ts_geom.xyz%', 'ts_final_geometry.xyz')
        replace_in_file(curr_orca_template, '%irc_F.xyz%', 'forward_end_opt.xyz')
        replace_in_file(curr_orca_template, '%irc_B.xyz%', 'backward_end_opt.xyz')
        replace_in_file(curr_orca_template, 'orca_Compound_1.hess', f'{basename}_Compound_1.hess')
        orca_basename = re.split(r'\.', os.bath.basename(curr_orca_template))[0]
        os.sy***REMOVED***m(f'orca.run {curr_orca_template} > {orca_basename}.out 2> {orca_basename}.err')
        os.chdir(init_path)
    os.chdir(init_path_top)
    

#######################################################################
#######################################################################
#######################################################################

if __name__== '__main__':
    print(' is module, not a script!')