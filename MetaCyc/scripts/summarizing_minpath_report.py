# Daniel Castaneda Mogollon, PhD
# February 4th, 2026
# Purpose:  This script will summarize the MinPath reports to print the number of Metabolic Pathways estimated by the
#           naive approach and by the parsimonious approach for each genome. It will also include the category of each pathway.
#           the summary of each pathway (if available), and the number of EC families found vs the expected number of EC families for each pathway.
# Input:    A report file(s) from MinPath, a detailed report file(s) from MinPath, and the Master MetaCyc database created previously with Category, Summary and other info.
# Output:   A summary table in tsv format with the information mentioned above for all genomes 

#%%
import pandas as pd
import os
import glob as glob

#%%
path = "/Users/danielcm/Desktop/SickKids/MinPath_2026/"
os.chdir(path)

master_df = pd.read_csv("Master_Metacyc_pathway_file.txt", sep = "\t")
for report in glob.glob("report/*.txt"):
    print(f"Processing file: {report}")
    
# %%
