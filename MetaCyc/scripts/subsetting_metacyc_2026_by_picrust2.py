# Daniel Castaneda Mogollon, PhD
# February 13th, 2026
# Purpose: This script was made to filter out the MetaCyc pathways from PICRUSt2
# so I can run MinPath with the database from PICRUSt2 and compare it to the default of MinPath
# and the latest MetaCyc 2026 database

#%%
import pandas as pd
import os

#%%
master_file_2026 = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/smart_tables_parsed/Master_Metacyc_pathway_file.tsv", sep="\t")

file_name = "/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_path2rxn_struc_filt_pro.txt"
picrust2_input = open(file_name, "r")
picrust2_db_formatted = open(file_name.replace(".txt","_formatted.txt"), "w")
picrust2_pwys = []
f_out = open("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/picrust2_pathways_reactions_subset.txt", "w")

pwy_and_reaction = {}

#icrust2_pwys = picrust2_input[0].tolist()
with picrust2_input as f:
    for line in f:
        line = line.strip()
        line = line.split("\t")
        pwy_id = line[0]
        picrust2_pwys.append(pwy_id)

for element in picrust2_pwys:
    if (master_file_2026.loc[:,"Pathways"] == element).any():
        pwy_and_reaction[element] = master_file_2026.loc[master_file_2026["Pathways"]==element,"Reactions of pathway"].values[0]
    
    else:
        print(f"Error: pathway {element} is in the PICRUSt2 database but not in the MetaCyc 2026 database. Exiting now.")
    #exit()

for item in pwy_and_reaction.keys():
    reaction = pwy_and_reaction[item].split(" // ")
    for element in reaction:
        f_out.write(item + "\t" + element + "\n")
        
print(f_out)
# %%
