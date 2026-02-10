# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:14:06 2022

@author: MSI
"""
from glob import glob
import pandas as pd
import re
import os
import numpy as np
from itertools import islice

if not os.path.exists('Human_metabolic_enzymes_AA_composition/'):
    os.mkdir('Human_metabolic_enzymes_AA_composition/')

folders = glob('AA_sequence/*')
for folder in folders:
    print('Processing of:  '+folder)
    files = glob(folder + '/*.txt')
            
    res = pd.DataFrame()
    ls_pro = []
    
    for f in files:
        pro = re.split(r'[/|\\|.]',f)[2]
        pathway = re.split(r'[/|\\|.]',f)[1]
        d = {}
        with open(f, 'r') as fin:
            lines = fin.readlines()
        if not lines:
            continue
        for line in lines[1:]:
            line = line.strip()
            for i in line:
                if i in d.keys():
                    d[i] += 1
                else:
                    d[i] = 1
                        
        AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
        
        #find the exceptional aa and delete it
        for i in list(d.keys()):
            if i not in AA:
                del d[i]
        
        d_sorted = dict(sorted(d.items(), key = lambda item:item[0]))
                    
        aa = list(d_sorted) 
        vacancy = list(set(AA) ^ set(aa))
        for i in vacancy:
            d_sorted[i] = 0
        
        d_resorted = dict(sorted(d_sorted.items(), key = lambda item:item[0])) 
        
        to = sum(list(d_sorted.values())) 
        
        for k in d_resorted.keys():
            d_resorted[k] = d_resorted[k]/to*100
            
        df = pd.DataFrame([d_resorted])
        res = res.append(df)
        ls_pro.append(pro)
        print(pro)
        
    res.index = ls_pro
    res.to_excel('Human_metabolic_enzymes_AA_composition/'+pathway+'.xlsx')
