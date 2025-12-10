# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 23:03:32 2024

@author: MSI
"""

import pandas as pd
from glob import glob
import numpy as np


files = glob('Human_metabolic_enzymes_AA_composition/*') 

df_results = pd.DataFrame()
for file in files:
    print('Processing of:  ' + file)
    df = pd.read_excel(file, index_col=(0))
    
    pathway_mean = df.mean(axis=0)
    pathway_mean = pd.DataFrame(pathway_mean).T
    df_results = df_results.append(pathway_mean)

df_results.to_excel('AA_pathways.xlsx', index=False)
