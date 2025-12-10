# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 10:22:48 2022

@author: MSI
"""
# In[0] Get pathway_list
import requests
import re
import pandas as pd

species = 'hsa'
link = 'https://rest.kegg.jp/list/pathway/' + species

req = requests.get(link)
req_text = req.text.split('\n')

cleaned_list = [item for item in req_text if item.strip()]

# Sort based on the number after 'hsa'
sorted_list = sorted(cleaned_list, key=lambda x: int(x.split('hsa')[1].split()[0]))

pathway_list = []
for item in sorted_list:
    pathway_list.append(item)
    if "hsa01250" in item:
        break
    
pathway_list = [item for item in pathway_list if "hsa01100" not in item] 

# In[0.5] Get full name of pathways
import os

if not os.path.exists('hsa/'): 
    os.mkdir('hsa/') 

pathway_name_list = [item.split('\t')[1] for item in pathway_list]
df_pathway_name = pd.DataFrame(pathway_name_list)
df_pathway_name.columns = ['Pathway']
df_pathway_name.to_excel('hsa/pathway_name.xlsx', index = False)


# In[1] Get pathway_link
import re
import pandas as pd
import numpy as np
import os

if not os.path.exists('hsa/'): 
    os.mkdir('hsa/') 

ls_pathway = []    
ls_pathway_link = []

for i in pathway_list:
    tmp = re.split(r'[:|\t]', i)
    pathway = tmp[0]
    print(pathway)
    ls_pathway.append(pathway)
    pathway_link = 'https://rest.kegg.jp/link/' + 'hsa/' + pathway 
    print(pathway_link)
    ls_pathway_link.append(pathway_link)

df = pd.DataFrame(ls_pathway_link)
df.to_excel('hsa/pathway_link.xlsx', index = False, header=None) 


# In[2] Make dir
import os

if not os.path.exists('hsa/AA_seq/'): 
    os.mkdir('hsa/AA_seq/') 
    
for p in ls_pathway:
    if not os.path.exists('hsa/AA_seq/'+p): 
            os.mkdir('hsa/AA_seq/'+p) 
            
# In[3] Crawling
import requests
import re

for link in ls_pathway_link:
    print('Processing of:  '+link)
    req = requests.get(link)
    req_txt = req.text.split('\n')
    req_txt = [item.replace("hsa:", "\thsa:") for item in req_txt] 
    while '' in req_txt:
        req_txt.remove('')
    for pp in req_txt:
        tmp1 = re.split(r'[\t]', pp)
        pathway1 = tmp1[0].split('path:')[1]
        protein = tmp1[2].split('hsa:')[1] 
        req_sequence = requests.get('https://rest.kegg.jp/get/'+'hsa:'+protein+'/aaseq')
        req_sequence_txt = req_sequence.text
        with open('hsa/AA_seq/'+pathway1+'/'+protein+'.txt','w') as f: 
            f.write(req_sequence_txt)
            print(pathway1+':'+protein)



