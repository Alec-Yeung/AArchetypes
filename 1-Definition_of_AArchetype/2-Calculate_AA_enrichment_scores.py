# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 19:26:06 2022

@author: MSI
"""

# In[1] Combine
from glob import glob
import pandas as pd

files = glob("Human_metabolic_enzymes_AA_composition/*")

combine = pd.DataFrame()
for xlsx in files:
    df = pd.read_excel(xlsx, index_col=(0))
    if len(df) < 10: 
        print(f"{xlsx} number of row is lower than 10，skip")
        continue
    combine = pd.concat([combine, df])
    
# In[2] Permutation test
from glob import glob
import pandas as pd
import numpy as np
import re

protein_pool = combine
protein_pool_index = protein_pool.index.values
protein_pool = protein_pool.reset_index()
n_pool = protein_pool.shape[0]

AA_pathways = pd.read_excel('AA_pathways.xlsx')
AA_pathways.drop(columns=['Pathway'], inplace=True)

d = {}
for xlsx in files:
    print('Processing of:  ' + xlsx)
    df = pd.read_excel(xlsx, index_col=(0))
    n_protein = len(df)
    
    if n_protein < 10: 
        print(f"{xlsx} number of row is lower than 10，skip")
        continue

    pathway = re.split(r'[\\|.]', xlsx)[1]
    real_AA_pathway = AA_pathways.loc[AA_pathways['Pathway_ID']==pathway].iloc[:,1:]
    real_AA_pathway = real_AA_pathway.values
    real_AA_pathway = real_AA_pathway.flatten()
    
    container = np.zeros((1000, 20))
    for i in range(1000):
        ran = np.random.random(n_protein)
        ran /= ran.max()
        index = np.array([int(r*n_pool) for r in ran])
        index[index >= n_pool-1] = n_pool-1
        virtual_pathway = protein_pool.iloc[index]

        virtual_pathway_values = virtual_pathway.iloc[:,1:].values
        virtual_pathway_mean = virtual_pathway_values.mean(axis=0)
        
        container[i] = virtual_pathway_mean
        
    p_list = []
    for u in range(20):
        x = sum(real_AA_pathway[u] > container[...,u])
        p = x/1000
        p_list.append(p)
    d[pathway] = p_list

df1 = pd.DataFrame(d).T
aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
df1.columns = aa

# In[3]
df2 = pd.read_excel('p_hsa.xlsx',index_col=(0)) 
df2[(df2>0.05) & (df2<0.95)] = 0.5  
df2.to_excel('p_adjust_hsa.xlsx') 
