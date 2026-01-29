# Daniel Castaneda Mogollon, PhD
# January 8th, 2026
# Script to convert MetaCyc pathway data to Anvio-compatible format by generating modules

import pandas as pd
import os

path = "/Users/danielcm/Desktop/diammatics/T1D/Anvio"
os.chdir(path)
metacyc_pathways = pd.read_csv("metacyc_rxn_to_level4ec.tsv", sep = "\t")
ec_file = pd.read_csv("ec_level4_info.tsv", sep = "\t")
#print(metacyc_pathways)

#I had to manually add the headers "Reaction" and "Enzyme 1 ... Enzyme 7" for pandas to read this file properly

module_dict = {}
description_dict = {}

for idx,row in metacyc_pathways.iterrows():
    reaction = row[0]
    module_dict[reaction] = []
    for enzyme in row[1:]:
        if pd.notna(enzyme):
            module_dict[reaction].append(enzyme)
            
for idx,row in ec_file.iterrows():
    ec_number = row[0]
    description = row[1]
    description_dict[ec_number] = description

not_found_counter = 0
for reaction in module_dict:
    module_name = "Metacyc_module_"+reaction
    with open(f"{module_name}.tsv", "w") as f:
        f.write(f"enzyme"+"\t"+"source"+"\t"+"orthology"+"\n")
        for enzyme in module_dict[reaction]:
            if enzyme not in description_dict:
                not_found_counter = not_found_counter + 1
                f.write(f"{enzyme}"+"\t"+"EGGNOG_EC_NUMBER"+"\t"+"Not in description file"+"\n")
            else:
                f.write(f"{enzyme}"+"\t"+"EGGNOG_EC_NUMBER"+"\t"+f"{description_dict[enzyme]}"+"\n")

print(f"Number of enzymes not found in description file: {not_found_counter}")
print("Modules generated successfully.")