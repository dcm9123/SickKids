import pandas as pd
from Bio import SeqIO
import os

path = "/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/MASTER/Mice/PICRUSt2.4/ns6_input_files/"
os.chdir(path)
comm = "ns6"
file = "NS6_filtered_ASVs_count600_len400_prev20.csv"


df = pd.read_csv(file, sep = ',')
sequence = df['asv_seq'].to_list()
f_out_name = file.replace(".txt",".fasta")
f_out = open(f_out_name,"w")
for seq in sequence:
    f_out.write(f">{seq}\n{seq}\n")

