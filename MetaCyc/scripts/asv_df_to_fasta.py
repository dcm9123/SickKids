# Daniel Castaneda Mogollon, PhD
# April 14th, 2026
# Purpose: This script takes as input the filtered ASV count tables, and takes the first. 
# column with the ASV IDs in the format of a sequence. Then it generates a fasta file with those
# IDs as the header and the sequence.


#%%
import pandas as pd
import re
from Bio import SeqIO
import glob as glob
import os

# %%
path = "/Users/danielcm/Desktop/Sickkids/Phyloseq2/Filtered_to_use/"
os.chdir(path)
comm = ["NS1","NS6","S2","S5"]
for file in glob.glob("*prev20*SRA*.csv"):
    df = pd.read_csv(file, sep = ',')
    sequence = df['asv_id'].to_list()
    f_out_name = file.replace(".csv",".fasta")
    f_out = open(f_out_name,"w")
    for seq in sequence:
        f_out.write(f">{seq}\n{seq}\n")

f_out.close()
# %%
