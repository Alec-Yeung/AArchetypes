# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 21:43:08 2024

@author: MSI
"""

# In[1] Entrez 2 Uniprot
import pandas as pd
import mygene

entrez = pd.read_excel('1100_Entrez_ID.xlsx') 

entrez_ls = entrez['Entrez_ID'].tolist()

mg = mygene.MyGeneInfo()

result = mg.querymany(entrez_ls, scopes='entrezgene', fields='uniprot', species='human')

uniprot_mapping = []
for entry in result:
    if 'uniprot' in entry:
        if isinstance(entry['uniprot'], dict) and 'Swiss-Prot' in entry['uniprot']:
            uniprot_mapping.append({'Entrez_ID': entry['query'], 'UniProt_ID': entry['uniprot']['Swiss-Prot']})
        elif isinstance(entry['uniprot'], list):
            for uniprot_entry in entry['uniprot']:
                if 'Swiss-Prot' in uniprot_entry:
                    uniprot_mapping.append({'Entrez_ID': entry['query'], 'UniProt_ID': uniprot_entry['Swiss-Prot']})

uniprot_df = pd.DataFrame(uniprot_mapping)

# In[2] Get PDB data
import os
import requests
from tqdm import tqdm

output_dir = "Alphafold_pdb_files"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

uniprot_ids = uniprot_df['UniProt_ID'].tolist()

base_url = "https://alphafold.ebi.ac.uk/files/"

for uniprot_id in tqdm(uniprot_ids):
    pdb_url = f"{base_url}AF-{uniprot_id}-F1-model_v4.pdb"
    
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        pdb_file_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
        with open(pdb_file_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded: {pdb_file_path}")
    else:
        print(f"Failed to download: {pdb_url}")

print("PDB files download complete!")

# In[3] Get the unobtained 
import os

obtained_ids = os.listdir("Alphafold_pdb_files/")
obtained_ids = [item.replace(".pdb", "") for item in obtained_ids]
obtained_ids = [item for item in obtained_ids if not isinstance(item, list)]

all_ids = uniprot_df['UniProt_ID'].tolist()
all_ids = [item for item in all_ids if not isinstance(item, list)]

unobtained_ids = list(set(all_ids) - set(obtained_ids)) 









