#!/usr/bin/env python
# coding: utf-8

# **STRUCTURE COMBINING**

# In[2]:


import os
f***REMOVED*** collections import defaultdict
f***REMOVED*** rdkit import Chem
f***REMOVED*** rdkit.Chem import RWMol, AllChem, rdMolDescriptors, Draw, rdmolfiles, rdDetermineBonds
f***REMOVED*** rdkit.Chem.rdmolfiles import MolsF***REMOVED***CDXMLFile, MolToXYZBlock, MolF***REMOVED***XYZFile


# In[7]:


mol = rdmolfiles.MolF***REMOVED***XYZFile('C15H12FNO2.xyz')
mol


# In[8]:


rdDetermineBonds.DetermineConnectivity(mol)
mol


# In[9]:


Chem.SanitizeMol(mol)
mol


# In[11]:


Chem.rdDetermineBonds.DetermineBondOrders(mol)
mol


# In[12]:


Chem.Kekulize(mol)
mol


# In[13]:


reference_smiles = 'FC=NC'
reference_fragment = Chem.MolF***REMOVED***Smiles(reference_smiles)
reference_fragment


# In[14]:


matches = mol.GetSubstructMatches(reference_fragment)
if not matches:
    raise ValueError("Reference fragment not found in the molecule!")
matches[0]  # Return the first match


# In[ ]:




