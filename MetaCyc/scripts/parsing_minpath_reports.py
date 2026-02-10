# Daniel Castaneda Mogollon, PhD
# February 5th, 2026
# Purpose: This script was made to parse each MinPath report and detailed report into one single
# file containing the info from all of the genomes analyzed.
# Input: The directories where the MinPath reports are located
# Output: A .tsv file with all the genomes and their predicted pathways by naïve approach, parsimonious, and the 
# default and updated MetaCyc databases.

#%%
import pandas as pd
import os
import glob as glob
# %%

dbs = ["default", "updated"]
dir1 = "/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/"
dir2 = "/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/detailed_reports/"


# %%

for method in dbs:
    pwy_dictionary = {}
    i = 0
    for file in glob.glob(os.path.join(dir1,method,"*.txt")):
        f_out = open(f"{file.split('/')[-1].split('_minpath')[0]}_{method}_parsed_report.tsv", "w")
        print(f"Processing file: {file}")
        with open(file,"r") as f:
            lines = f.readlines()
            pwy = line.split(" ")[1]
            for line in lines:
                i = i + 1
                line = line.strip()
                if "naive 1" in line:
                    pwy_dictionary[pwy] = {"naive":"Yes"}
                elif "naive 0" in line:
                    pwy_dictionary[pwy] = {"naive":"No"}
                if "minpath 1" in line:
                    pwy_dictionary[pwy] = {"minpath":"Yes"}
                elif "minpath 0" in line:
                    pwy_dictionary[pwy] = {"minpath":"No"}
                pwy_dictionary[pwy] = {"Total families": line.split("  ")[4]}
                pwy_dictionary[pwy] = {"Families found": line.split("  ")[6]}
            sum = len(naive) + len(not_naive) + len(minpath) + len(not_minpath)
            if i!=sum:
                exit(f"Error in file {file}: the sum of the categories does not match the total number of lines processed. Exiting now")
            for element,element2 in zip(fam,fam_found):
                print(element,element2)
                try:
                    int(element)
                    int(element2)
                except ValueError:
                    exit(f"Error in file {file}: the family counts are not integers. Exiting now")
        
            
        break  
# %%
