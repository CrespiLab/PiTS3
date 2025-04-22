#!/usr/bin/env python3

import os, re, shutil, sys, glob, math, subprocess
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdDetermineBonds, rdMolTransforms
from ts_pipeline.config import TOOLS as tsp_tools

### Installation section
pysis_path  = tsp_tools["pysis"]
xtb_path    = tsp_tools["xtb"]
crest_path  = tsp_tools["crest"]
orca_path   = tsp_tools["orca"]
ts_pipe_dir = os.path.dirname(os.path.realpath(__file__))

#%%
### CATEGORY: Fragment identification functions
def find_fragment_atoms(xyzfile, reference_smiles, chrg=0, sanitize=False):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolFromMol2Block(mol.write('mol2'), sanitize=sanitize)
    reference_fragment = Chem.MolFromSmiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    if len(matches) > 1:
        raise ValueError("Selected reference fragment ({reference_smiles}) is not unique in the molecule!")
    return matches[0]  # Return the first match    
def find_fragment_atoms_with_neighbors(xyzfile, reference_smiles, chrg=0, sanitize=False):
    """
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolFromMol2Block(mol.write('mol2'), removeHs = False, sanitize=sanitize)
    reference_fragment = Chem.MolFromSmiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    if len(matches) > 1:
        raise ValueError("Selected reference fragment (--dihedral) is not unique in the molecule!")
    key_atoms = [mol.GetAtomWithIdx(atom_index) for atom_index in matches[0][1:3]]
    neighbors_indices = []
    for atom in key_atoms:
        neighbors_idxs = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in matches[0]]
        neighbors_indices.extend(neighbors_idxs)
    all_atoms = neighbors_indices + list(matches[0])
    return all_atoms    
def get_dihedral(xyzfile, atom1, atom2, atom3, atom4):
    '''
    Measures dihedral angle in an .xyz molecule.
    
    Parameters
    ----------
    xyzfile : str
        xyz geometry file.
    atom1, atom2, atom3, atom4 : ints
        "Natural" atom numbers (with numbering starting at 1!)

    Returns
    -------
    Dihedral in degrees.

    '''
    try:
        mol = rdmolfiles.MolFromXYZFile(xyzfile)
        atoms = list(map(lambda x: x-1, [atom1, atom2, atom3, atom4]))
        dihedral = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *atoms)
        return dihedral
    except AttributeError as e:
        print(e)
        return None
def get_angle(xyzfile, atom1, atom2, atom3):
    '''
    Measures angle in an .xyz molecule.
    
    Parameters
    ----------
    xyzfile : str
        xyz geometry file.
    atom1, atom2, atom3 : ints
        "Natural" atom numbers (with numbering starting at 1!)

    Returns
    -------
    Angle in degrees.

    '''
    try:
        mol = rdmolfiles.MolFromXYZFile(xyzfile)
        atoms = list(map(lambda x: x-1, [atom1, atom2, atom3]))
        angle = rdMolTransforms.GetAngleDeg(mol.GetConformer(), *atoms)
        return angle
    except AttributeError as e:
        print(e)
        return None
def get_distance(xyzfile, atom1, atom2):
    '''
    Measures distance in an .xyz molecule.
    
    Parameters
    ----------
    xyzfile : str
        xyz geometry file.
    atom1, atom2 : ints
        "Natural" atom numbers (with numbering starting at 1!)

    Returns
    -------
    Distance in ångström (hello from Uppsala University!).

    '''
    try:
        mol = rdmolfiles.MolFromXYZFile(xyzfile)
        atoms = list(map(lambda x: x-1, [atom1, atom2]))
        distance = rdMolTransforms.GetBondLength(mol.GetConformer(), *atoms)
        return distance
    except AttributeError as e:
        print(e)
        return None
def _find_fragment_atoms_rdkit(xyzfile, reference_smiles, chrg=0):
    """
    Deprecated.
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.

    Fixes:
    - Uses RDKit's DetermineConnectivity() to assign bonds to XYZ molecules.
    - Sanitizes the molecule to ensure valence is properly assigned.
    
    :param mol: RDKit molecule
    :param reference_fragment: RDKit molecule representing the common fragment
    :return: List of atom indices that match the fragment
    """
    mol = rdmolfiles.MolFromXYZFile(xyzfile, )
    rdDetermineBonds.DetermineConnectivity(mol, charge=chrg)
    rdDetermineBonds.DetermineBondOrders(mol, charge=chrg)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    reference_fragment = Chem.MolFromSmiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    return matches[0]  # Return the first match    
def _find_fragment_atoms_with_hydrogens(xyzfile, reference_smiles, chrg=0, sanitize=False):
    """
    Deprecated.
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.
    Uses openbabel to generate mol object with correct bond orders.

    """
    mol = next(pybel.readfile('xyz', xyzfile))
    mol = rdmolfiles.MolFromMol2Block(mol.write('mol2'), removeHs = False, sanitize=sanitize)
    reference_fragment = Chem.MolFromSmiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    if len(matches) > 1:
        raise ValueError("Selected reference fragment (--dihedral) is not unique in the molecule!")
    key_atoms = [mol.GetAtomWithIdx(atom_index) for atom_index in matches[0][1:3]]
    hydrogen_indices = []
    for atom in key_atoms:
        hydrogen_idxs = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"]
        hydrogen_indices.extend(hydrogen_idxs)
    all_atoms = hydrogen_indices + list(matches[0])
    return all_atoms    
def _find_fragment_atoms_ibo(xyzfile, reference_smiles, chrg=0):
    """
    Deprecated.
    Finds and extracts the atom indices of a reference fragment in a molecule from an XYZ file.
    IGNORES BOND ORDERS

    Fixes:
    - Uses RDKit's DetermineConnectivity() to assign bonds to XYZ molecules.
    - Sanitizes the molecule to ensure valence is properly assigned.
    
    :param mol: RDKit molecule
    :param reference_fragment: RDKit molecule representing the common fragment
    :return: List of atom indices that match the fragment
    """
    mol = rdmolfiles.MolFromXYZFile(xyzfile, )
    rdDetermineBonds.DetermineConnectivity(mol, charge=chrg)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    reference_fragment = Chem.MolFromSmiles(reference_smiles)
    matches = mol.GetSubstructMatches(reference_fragment)
    if not matches:
        raise ValueError("Reference fragment not found in the molecule!")
    return matches[0]  # Return the first match    
#%%
### CATEGORY: xTB invoking functions
def xtb_opt(input_file, dirname='.', model='--gfn2', solvent='', optlev='', etemp=''):
    '''
    Optimizes molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - directory to run xTB job in
        Returns:
    opt_path - XYZ geometry of optimized molecule
    '''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    os.system(f'{xtb_path} {input_file} --opt {model} {solvent} {optlev} {etemp} | tee xtb_opt_stdout.log')
    opt_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return opt_path
def xtb_scan_rotation(input_file, dirname='.', model='--gfn2', solvent='',
                      dihedral='0,0,0,0,0', scan='0,0,0', optlev='', etemp=''):
    '''
    Run a scan along selected dihedral of molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    rotated_path - XYZ geometry of alternative diastereomer
    '''
    scan_input=f'''$constrain
 force constant=0.50
 dihedral: {dihedral}
$scan
 1: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.system(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} {etemp} | tee xtb_dih_scan_stdout.log')
    rotated_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return rotated_path
def xtb_scan_nbd(input_file, dirname='.', model='--gfn2', solvent='',
                      distance1='0,0', distance2='0,0', scan='0,0,0', optlev='', etemp=''):
    '''
    Run a concerted scan along distance between double bond in double bonds in norbornadiene at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    final_path - XYZ geometry of alternative diastereomer
    '''
    scan_input=f'''$constrain
   force constant=0.5
   distance: {distance1[0]}, {distance1[1]}, auto
   distance: {distance2[0]}, {distance2[1]}, auto
$scan
   mode=concerted
   1: {scan}
   2: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.system(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} {etemp} | tee xtb_dih_scan_stdout.log')
    final_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return final_path
def xtb_scan_adae(input_file, dirname='.', model='--gfn2', solvent='',
                      distance1='0,0', scan='0,0,0', optlev='', etemp=''):
    '''
    Run a concerted scan along distance between double bond in double bonds in norbornadiene at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    final_path - XYZ geometry of alternative diastereomer
    '''
    scan_input=f'''$constrain
   force constant=0.5
   distance: {distance1[0]}, {distance1[1]}, auto
$scan
   1: {scan}
$end'''
    input_file, initial_path = safe_dir(input_file, dirname, rename='initial_structure.xyz')
    with open('scan.inp','w') as file:
        file.write(scan_input)
    os.system(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} {etemp} | tee xtb_dih_scan_stdout.log')
    final_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return final_path
def _xtb_scan_inversion(input_file, dirname='.', model='--gfn2', solvent='',
                       angle='0,0,0,0.0', scan='0,0,0', optlev=''):
    '''
    Deprecated.
    INVERSION FOR MORE THAN 180 DEGREES IS PROBLEMATIC
    Run a scan along selected angle of molecule filename.xyz with xtb at gfn2 level of theory
        Arguments:
    path         - path to XYZ geometry file
        Returns:
    inverted_path - XYZ geometry of alternative diastereomer
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
    os.system(f'{xtb_path} {input_file} --opt --input scan.inp {model} {solvent} {optlev} | tee xtb_angle_scan.log')
    inverted_path = os.path.abspath('xtbopt.xyz')
    os.chdir(initial_path)
    return inverted_path
#%%
### CATEGORY: pysisyphus invoking functions
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
        sys.exit("You need exactly two geometries to continue, exiting")
    input_yaml, initial_path = safe_dir(input_yaml, dirname)    
    xyz_1 = safe_dir(xyz_1, dirname, rename='geom_1.xyz')[0]
    xyz_2 = safe_dir(xyz_2, dirname, rename='geom_2.xyz')[0]
    os.chdir(dirname)
    os.system(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
    shutil.rmtree('qm_calcs')
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
        os.system(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
        shutil.rmtree('qm_calcs')
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
        os.system(f'{pysis_path} {input_yaml} | tee pysis_stdout.log')
        shutil.rmtree('qm_calcs')
        reoptimized_TS_conformer = re.split('\.', curr_TS_conformer)[0] + '_reopt.xyz'
        # breakpoint()
        reoptimized_TSes['conformers'][f'{reoptimized_TS_conformer}'] = {'backward' : os.path.abspath('backward_end_opt.xyz'),
                                                               'TS' : os.path.abspath('ts_final_geometry.xyz'),
                                                               'forward' : os.path.abspath('forward_end_opt.xyz')}
        os.chdir(init_path)
    os.chdir(init_path_top)
    return reoptimized_TSes
#%%
### CATEGORY: CREST invoking functions    
def crest_constrained_sampling(input_file,
                               dirname='.',
                               model='--gfn2', 
                               solvent='', 
                               constrain_atoms=[],
                               constrain_dihedral=[],
                               force_constant = 0.50,
                               optlev='', 
                               dlen='',
                               mdlen='',
                               cinp='constraints.inp'):
    '''
    Runs CREST conformational sampling with constrain (predefined in the function).
    '''
    input_file, initial_path = safe_dir(input_file, dirname)    
    os.chdir(dirname)
    constrain_atoms = list(map(lambda x: x+1, constrain_atoms))
    constrain_atoms_line = ','.join(map(str,constrain_atoms))
    if constrain_dihedral:
        constrain_dihedral_line = ','.join(map(str,constrain_dihedral)) + ',auto'
        constraint_input=f'''$constrain
force constant={force_constant}
 atoms: {constrain_atoms_line}
 dihedral: {constrain_dihedral_line}
$end'''
    else:
        constraint_input=f'''$constrain
 force constant={force_constant}
 atoms: {constrain_atoms_line}
$end'''
    with open('constraints.inp','w') as file:
        file.write(constraint_input)
    crest_line=f'{crest_path} {input_file} --cinp {cinp} {optlev} {model} {solvent} {dlen} {mdlen} | tee xtb_TS_conf_sampling_stdout.log'
    os.system(crest_line)
    TS_conformers_path = f'{dirname}/crest_conformers.xyz'
    os.chdir(initial_path)
    return TS_conformers_path
def crest_simple_sampling(input_file,
                               dirname='.',
                               model='--gfn2', 
                               solvent='', 
                               optlev='', 
                               dlen='',
                               mdlen='',):
    '''
    Runs CREST conformational sampling with constrain (predefined in the function).
    '''
    input_file, initial_path = safe_dir(input_file, dirname)    
    os.chdir(dirname)
    crest_line=f'{crest_path} {input_file} {optlev} {model} {solvent} {dlen} {mdlen} | tee xtb_TS_conf_sampling_stdout.log'
    os.system(crest_line)
    conformers_path = f'{dirname}/crest_conformers.xyz'
    os.chdir(initial_path)
    return conformers_path
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
    os.system(f'{crest_path} {ensemble_file} --cregen {ensemble_file} {extra_params} | tee cregen_stdout.log')
    basename = os.path.basename(ensemble_file)
    basename = re.split(r'\.', basename)[0]
    if sorted_ensemble == True:
        resulting_ensemble = os.path.abspath(f'{basename}.xyz.sorted')
    else:
        resulting_ensemble = os.path.abspath('crest_ensemble.xyz')
    os.chdir(init_path)
    return resulting_ensemble
def _crest_constrained_sampling_old(input_file,
                                  dirname='.',
                                  model='--gfn2', 
                                  solvent='', 
                                  angles=[],
                                  dihedral=[],
                                  optlev='', 
                                  dlen='',
                                  mdlen='',
                                  cinp='constraints.inp'):
    '''
    Deprecated.
    Runs CREST conformational sampling with constrain (predefined in the function).
    '''
    input_file, initial_path = safe_dir(input_file, dirname)    
    os.chdir(dirname)
    dihedral = ','.join(map(str,dihedral))
    dihedral = dihedral + ', auto'
    constraint_input=f'''$constrain
 force constant=5.00
 dihedral: {dihedral}
 angle: {angles[0]}
$end'''
    if len(angles) > 1:
        constraint_input=f'''$constrain
 force constant=5.00
 dihedral: {dihedral}
 angle: {angles[0]}
 angle: {angles[1]}
    $end'''
    with open('constraints.inp','w') as file:
        file.write(constraint_input)
    crest_line=f'{crest_path} {input_file} --cinp {cinp} {optlev} {model} {solvent} {dlen} {mdlen} | tee xtb_TS_conf_sampling_stdout.log'
    os.system(crest_line)
    TS_conformers_path = f'{dirname}/crest_conformers.xyz'
    os.chdir(initial_path)
    return TS_conformers_path
#%%
### CATEGORY: ORCA invoking functions
def orca_three_points(irc_dict, orca_template = 'orca_three_points.inp', dirname = '.',
                      postpone_orca = False,):
    '''
    Runs ORCA compound job to reoptimize TS and IRC endpoints. DOES NOT DO DFT IRC!

    Parameters
    ----------
    irc_dict : str
        Dictionary from pysis_ts_irc, containing paths to TS and both endpoins for each TS conformer.
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
        orca_basename = re.split(r'\.', os.path.basename(curr_orca_template))[0]
        if not postpone_orca:
            os.system(f'{orca_path} {curr_orca_template} > {orca_basename}.out 2> {orca_basename}.err')
            for f in glob.glob('cregened*tmp*'):
                os.remove(f)
        else:
            print(f'--postpone_orca requested. ORCA input files and geometries have been prepared,\nbut you will have to start orca jobs manually.\nExiting now.')
        os.chdir(init_path)
    os.chdir(init_path_top)
def orca_multixyz(multixyzfile, 
                  orca_template = 'orca_single_point.inp', 
                  dirname = '.',
                  postpone_orca = False,):
    '''
    Runs ORCA compound job to reoptimize TS and IRC endpoints. DOES NOT DO DFT IRC!

    Parameters
    ----------
    irc_dict : str
        Dictionary from pysis_ts_irc, containing paths to TS and both endpoins for each TS conformer.
    orca_template : str
        ORCA compound job template. Do not change if you do not know how to set compound job!
    dirname : TYPE, optional
        Directory in which to run the job. The default is '.'.

    Returns
    -------
    Path to orca output.

    '''
    init_path_top = os.getcwd()
    multixyzfile = os.path.abspath(multixyzfile)
    orca_template = os.path.abspath(orca_template)
    os.chdir(dirname)
    multixyz = shutil.copy(multixyzfile, '.')
    multixyz = split_geoms(multixyz)
    for file in multixyz[0]:
        name = re.split(r'\.', file)[0]
        conf_dir = name + '_orca_SP'
        os.mkdir(conf_dir)
        curr_file = shutil.move(file, conf_dir)
        curr_orca_template = shutil.copy(orca_template, conf_dir)
        replace_in_file(curr_orca_template, '%geom.xyz%', f'{os.path.basename(curr_file)}')
        orca_basename = re.split(r'\.', os.path.basename(curr_orca_template))[0]
        if not postpone_orca:
            subprocess.run(f'{orca_path} {os.path.basename(curr_orca_template)} > {orca_basename}.out 2> {orca_basename}.err', cwd = f'{conf_dir}',
                           shell=True, capture_output = True, text = True)
            for f in glob.glob('cregened*tmp*'):
                os.remove(f)
        else:
            print(f'--postpone_orca requested. ORCA input files and geometries have been prepared,\nbut you will have to start orca jobs manually.\nExiting now.')
    os.chdir(init_path_top)
def _orca_three_points_with_control(irc_dict, orca_template = 'orca_three_points.inp', dirname = '.',
                      postpone_orca = False,
                      control_ang=[], control_ang_range=[],
                      control_dih=[], control_dih_range=[]):
    '''
    TODO
    Runs ORCA compound job to reoptimize TS and IRC endpoints. DOES NOT DO DFT IRC!

    Parameters
    ----------
    irc_dict : str
        Dictionary from pysis_ts_irc, containing paths to TS and both endpoins for each TS conformer.
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
        if control_ang:
            control_ang_value = get_angle(xyzs['TS'], *control_ang)
            control_ang_range.sort()
            if control_ang_range[0] < control_ang_value < control_ang_range[1]:
                print(f'Angle {control_ang} is in undesired range {control_ang_range}, skipping {conformer}')
                continue
        if control_dih:
            control_dih_value = get_dihedral(xyzs['TS'], *control_dih)
            control_dih_range.sort()
            if control_dih_range[0] < control_dih_value < control_dih_range[1]:
                print(f'Angle {control_dih} is in undesired range {control_dih_range}, skipping {conformer}')
                continue
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
        orca_basename = re.split(r'\.', os.path.basename(curr_orca_template))[0]
        if not postpone_orca:
            os.system(f'{orca_path} {curr_orca_template} > {orca_basename}.out 2> {orca_basename}.err')
            for f in glob.glob('cregened*tmp*'):
                os.remove(f)
        else:
            print(f'--postpone_orca requested. ORCA input files and geometries have been prepared,\nbut you will have to start orca jobs manually.\nExiting now.')
        os.chdir(init_path)
    os.chdir(init_path_top)
#%%
### CATEGORY: auxiliary functions
def mkbasedir(path_to_xyzfile, prefix='', suffix='', ignore_existing=False):
    '''
    Creates a directory with the basename of a file
        Arguments:
    filename (str) - name of the file to take basename from
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

def _step():
    '''
    Despecated
    Automatic step counter
    '''
    def number_generator():
        num = 1
        while True:
            yield num
            num += 1
    gen = number_generator()  # Create the generator instance
    def next_number():
        print('Step ', next(gen), ' of TS_pipe script')  # Automatically print the next number
    return next_number  # Return the function that prints the next number
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
def orca_user_confirmation():
    while True:
        user_input = input("No ORCA template provided for this run, thus it will stop after IRC is done for all TS conformers. Are you sure? (yes/no) [yes]: ")
        if user_input == "":
            user_input = "yes"
        if user_input in ("yes", "y", "no", "n"):
            break
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")
    if user_input in ("no", "n"):
        sys.exit('Stopping now because no ORCA template provided, and user requested stop')
def confirm_orca_template_exists(orca_template_file, ask_yes = True):
    while True:
        if os.path.isfile(orca_template_file):
            break
        else:
            print(f'ORCA template {orca_template_file} not found, trying to find similar template in ts_pipeline/templates...')
            template_dir = os.path.abspath(f'{ts_pipe_dir}/../templates')
            if os.path.isfile(f'{template_dir}/{os.path.basename(orca_template_file)}') and ask_yes:
                confirmation = input(f'Found {os.path.basename(orca_template_file)} in {template_dir}, confirm? [yes] ')
                if 'y' in confirmation or not confirmation:
                    return os.path.abspath(f'{template_dir}/{os.path.basename(orca_template_file)}')
                elif 'n' in confirmation:
                    sys.exit('Aborted by user (ORCA template not found)')
            elif os.path.isfile(f'{template_dir}/{os.path.basename(orca_template_file)}') and not ask_yes:
                print(f'{orca_template_file} not found, {os.path.basename(orca_template_file)} found instead in {template_dir}. Taking the last one as template.')
                return os.path.abspath(f'{template_dir}/{os.path.basename(orca_template_file)}')
            elif not os.path.isfile(f'{template_dir}/{os.path.basename(orca_template_file)}'):
                sys.exit(f'ORCA template suggested by user, but no template file found in current directory nor among templates. Aborting.')
#%%
### CATEGORY: processing functions                    
def Williams_proc1(mol, e_dict, verbose = False):
#    R = 8.31446261815324
    R = 8.314
    T = 298.15
    hartree = 2625.4996394799 # [kJ/mol]
    min_energy = min(min(e_dict["Metastable energies"]), min(e_dict["TS energies"]))
    M_Es = [(M_E - min_energy) * hartree for M_E in e_dict['Metastable energies']]
    TS_Es = [(TS_E - min_energy) * hartree for TS_E in e_dict['TS energies']]
    M_BWs = [math.exp(-1000*G/ (R*T)) for G in M_Es]
    TS_BWs = [math.exp(-1000*G/ (R*T)) for G in TS_Es]
    TS_BWs_sum = sum(TS_BWs)
    khi_zero = max(M_BWs) / sum(M_BWs)
    BWeff = khi_zero * TS_BWs_sum
    Geff = (-R * T * math.log(BWeff)) / 1000
    if verbose:
        print(f'\nWilliams (10.1002/poc.4312) procedure 1 for {mol}\n'+f'-------------------|'*3)
        print(f'Energies [hartree] | Energies [kJ/mol] | Boltzmann weights |')
        print(f'   Mi        TSi   |'*3)
        print(*[f'{M_Eh: <10.4f}{TS_Eh: <9.4f}|{M_E: ^10.1f}{TS_E: ^9.1f}|{M_BW: <9.3e}|{TS_BW: <9.3e}|' 
              for M_Eh, TS_Eh, M_E, TS_E, M_BW, TS_BW 
              in zip(e_dict['Metastable energies'], e_dict["TS energies"], M_Es, TS_Es, M_BWs, TS_BWs)], sep='\n')
        print(f'-------------------|'*3)
        print(f'Q(TS): {TS_BWs_sum:.3e}, χ0 = {khi_zero:.3e}, χ0 * Q(TS) = {TS_BWs_sum * khi_zero:.3e}')
        print(f'Geff({mol}) = {Geff:.3f} kJ/mol')
    return Geff
def Williams_proc3(mol, e_dict, verbose = False):
#    R = 8.31446261815324
    R = 8.314
    T = 298.15
    hartree = 2625.4996394799 # [kJ/mol]
    min_energy = min(min(e_dict["Metastable energies"]), min(e_dict["TS energies"]))
    M_Gs = [(M_G - min_energy) * hartree for M_G in e_dict['Metastable energies']]
    TS_Gs = [(TS_G - min_energy) * hartree for TS_G in e_dict['TS energies']]
    M_BWs = [math.exp(-1000*G/ (R*T)) for G in M_Gs]
    TS_BWs = [math.exp(-1000*G/ (R*T)) for G in TS_Gs]
    Q_TS = sum(TS_BWs)
    Q_M = sum(M_BWs)
    Geff = R * T * math.log(Q_M / Q_TS) / 1000
    if verbose:
        print(f'\nWilliams (10.1002/poc.4312) procedure 3 for {mol}\n'+f'-------------------|'*3)
        print(f'   Gibbs [hartree] | Rel.Gibbs[kJ/mol] | Boltzmann weights |')
        print(f'   Mi        TSi   |'*3)
        print(*[f'{M_Gh: <10.4f}{TS_Gh: <9.4f}|{M_G: ^10.1f}{TS_G: ^9.1f}|{M_BW: <9.3e}|{TS_BW: <9.3e}|' 
              for M_Gh, TS_Gh, M_G, TS_G, M_BW, TS_BW 
              in zip(e_dict['Metastable energies'], e_dict["TS energies"], M_Gs, TS_Gs, M_BWs, TS_BWs)], sep='\n')
        print(f'-------------------|'*3)
        Gavg_M = R * T * sum([(math.exp(-1000 * G / ( R * T )) / Q_M ) * 1000 * G / ( R * T ) for G in M_Gs])
        Smix_M = R * T * sum([ (math.exp(-1000 * G / (R * T)) / Q_M ) * math.log(math.exp(-1000 * G / (R * T)) / Q_M) for G in M_Gs])
        Geff_M = Gavg_M + Smix_M
        print(f'Q(M)  = {Q_M: .3e}, Gavg(M)  = {Gavg_M: .3e}, ΔSmix(M)  = {Smix_M: .3e}, Geff(M)  = {Geff_M: .3e}', sep='')
        Gavg_TS = R * T * sum([(math.exp(-1000 * G / ( R * T )) / Q_TS ) * 1000 * G / ( R * T ) for G in TS_Gs])
        Smix_TS = R * T * sum([ (math.exp(-1000 * G / (R * T)) / Q_TS ) * math.log(math.exp(-1000 * G / (R * T)) / Q_TS) for G in TS_Gs])
        Geff_TS = Gavg_TS + Smix_TS
        print(f'Q(TS) = {Q_TS: .3e}, Gavg(TS) = {Gavg_TS: .3e}, ΔSmix(TS) = {Smix_TS: .3e}, Geff(TS) = {Geff_TS: .3e}', sep='')
        print(f'Geff({mol}) = {Geff:.3f} kJ/mol')
    return Geff
#######################################################################
#######################################################################
#######################################################################

if __name__== '__main__':
    print('This is module, not a script!')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    