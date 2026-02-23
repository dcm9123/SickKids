# Daniel Castaneda Mogollon, PhD
# February 11th, 2026
# Purpose: This script was made to filter out the MetaCyc pathways from PICRUSt2 so I can annotate them myself using
# MinPath. The annotation will be made on the genomes.

#%%
import pandas as pd
import os

#%%
picrust2_db = pd.read_csv("/Users/danielcm/anaconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_pro.txt", sep = "\t", header=None)
metacyc_2026 = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/smart_tables_parsed/Master_Metacyc_pathway_file.tsv", sep="\t")
pathways_metacyc_2026 = list(metacyc_2026['Pathways'])
picrust2_pathways = list(picrust2_db[0])

overlapping_pwys = []
overlapping_ecs = []

for row in metacyc_2026.itertuples():
    if row[1] in picrust2_pathways:
        overlapping_pwys.append(row[2])
        overlapping_ecs.append(row[2])

for pwy in picrust2_pathways:
    if pwy not in pathways_metacyc_2026:
        print(f"Error: pathway {pwy} is in the overlapping pathways list but not in the PICRUSt2 database. Exiting now.")
        #exit()
        
print(f"The PICRUSt2 database has a total of: {len(picrust2_pathways)} pathways")
print(f"The MetaCyc 2026 database has a total of: {len(pathways_metacyc_2026)} pathways")
print(f"The number of pathways that overlap between the two databases is: {len(overlapping_pwys)}")
print(f"The number of pathways that overlap between the two databases is: {len(set(overlapping_pwys))}")

#print(overlapping_pwys)

# %%
