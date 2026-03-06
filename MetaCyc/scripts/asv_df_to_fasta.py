#%%
import pandas as pd
import re
from Bio import SeqIO
import glob as glob

# %%
path = "/Users/danielcm/Desktop/Sickkids/Phyloseq2/"
os.chdir(path)
comm = ["ns1","ns6","s2","s5"]
for file in glob.glob("*0.txt"):
    df = pd.read_csv(file, sep = '\t')
    sequence = df['ASV_ID'].to_list()
    f_out_name = file.replace(".txt",".fasta")
    f_out = open(f_out_name,"w")
    for seq in sequence:
        f_out.write(f">{seq}\n{seq}\n")

f_out.close()
# %%
