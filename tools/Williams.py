# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 19:03:38 2025

@author: ***REMOVED***st136
"""

# import os, sys, re, argparse, pdb, glob, pprint, json
# import TS_pipe as tsp
import pandas as pd
import numpy as np
import scipy.constants as Constants
# f***REMOVED*** rdkit import Chem
# f***REMOVED*** rdkit.Chem import rdmolfiles
# f***REMOVED*** openbabel import pybel

R = Constants.R
T = 298.15 # default temperature K
Constant_Hartree_to_kJmol = 2625.4996394799

def Hartree_to_kJmol(Hartrees):
    result = Hartrees * Constant_Hartree_to_kJmol
    return result

def exp_RT(deltaG):
    '''
    Parameters
    ----------
    deltaG : TYPE
        DESCRIPTION.

    Returns
    -------
    result : TYPE
        DESCRIPTION.

    '''
    result = np.exp(-deltaG/(R*T))
    return result

def expRT_to_deltaG(expRT_Geff):
    '''
    Parameters
    ----------
    expRT_Geff : TYPE
        DESCRIPTION.

    Returns
    -------
    result : TYPE
        DESCRIPTION.

    '''
    result = -R*T*np.log(expRT_Geff)
    return result

#%%


# json_file = './stilbenes_final_energies.json'
# df = pd.read_json(json_file, orient='index')
# df = df.reset_index(names='Molecule')
# df



#%% TEST

df = pd.DataFrame({"Z energies":[[-978.760013,-978.758859,-978.758977,-978.758864]],
                   "TS energies": [[-978.724000,-978.724958,-978.723867,-978.724511]]
                   })
df["Z energies"] = df["Z energies"].apply(lambda x: Hartree_to_kJmol(np.array(x)))
df["TS energies"] = df["TS energies"].apply(lambda x: Hartree_to_kJmol(np.array(x)))
df

#%% CALCULATE SUM of Gibbs energies of activation
df['Z energies (refd)'] = df['Z energies'].apply(lambda x: np.array(x) - np.min(x)) # subtract lowest Z energy f***REMOVED*** all
df['Z energy lowest'] = df['Z energies'].apply(min) # get lowest Z energy

df["TS minus Zlowest"] = df.apply(lambda row: np.array(row["TS energies"]) - row["Z energy lowest"], axis=1) # subtract lowest Z f***REMOVED*** all TS
# df["TS minus Zlowest"]

df["TS minus Zlowest Jmol"] = df["TS minus Zlowest"]*1000


df["exp_RT of TS-Zlowest"] = df['TS minus Zlowest Jmol'].apply(lambda x: exp_RT(x)) # calculate exp_RT for all TS-Zlowest
# df["exp_RT of TS-Zlowest"]
df["SUM exp_RT TS-Zlowest"] = df["exp_RT of TS-Zlowest"].apply(sum) # calculate sum of exp_RT
# df["SUM exp_RT"]

#%% CALCULATE MOLE FRACTION x0
df['exp_RT of Zrefd'] = df['Z energies (refd)'].apply(lambda x: exp_RT(x))
# df['exp_RT of Zrefd']
df['SUM exp_RT Zrefd'] = df['exp_RT of Zrefd'].apply(sum)
# df['SUM exp_RT Zrefd']

df['exp_RT Zlowest'] = df['Z energies (refd)'].apply(min).apply(lambda x: exp_RT(x)) # exp_RT(0) = 1
# df['exp_RT Zlowest']
df['x_0_RS'] = df.apply(lambda row: row['exp_RT Zlowest'] / row['SUM exp_RT Zrefd'], axis=1) # calculate mole fraction of x_0
# df['x_0_RS']

df['expRT Geff'] = df.apply(lambda row: row['x_0_RS'] * row['SUM exp_RT TS-Zlowest'], axis=1)
df['expRT Geff']

df['deltaGact'] = df['expRT Geff'].apply(lambda x: expRT_to_deltaG(x))
df['deltaGact']

df['deltaGact_kJmol'] = df['deltaGact']/1000
df['deltaGact_kJmol']


#%%