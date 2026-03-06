# Daniel Castaneda Mogollon, PhD
# February 21st, 2026
# Purpose: To crossmap and convert each ASV from PICRUSt2.3 to its corresponding taxa from PICRUSt2.3 run

#%%
import pandas as pd
import os
import glob as glob
from Bio import SeqIO
from Bio.SeqIO import FastaIO
#%%

path = "/Users/danielcm/Desktop/SickKids/"
os.chdir(path)

asv_dict = {}
communities = ["ns1","ns6","s2","s5"]

df1 = pd.read_csv("Phyloseq2/samples_364_by_seq_taxa_table_genus_species.csv",sep = ",")
list_id = []
set_id = set()
list_seq = []
set_seq = set()

'''
for com in communities:
    community_fasta = open(f"Phyloseq/ps_{com}_asv_sequences.fasta", 'r')
    for record in SeqIO.parse(community_fasta, "fasta"):
        asv_id = record.id
        asv_seq = str(record.seq)
        asv_dict[asv_id] = asv_seq   
        list_id.append(asv_id)
        set_id.add(asv_id)
        list_seq.append(asv_seq)
        set_seq.add(asv_seq)

print(f"There are a total of {len(asv_dict.keys())} unique ASV IDs across the four consortia, confirmed by a Python dictionary.")
print(f"There are a total of {len(set_id)} unique ASV IDs across the four consortia, confirmed by a Python set.")
print(f"There are a total of {len(set_seq)} unique ASV sequences across the four consortia, confirmed by a Python set.")
print(f"There are a total of {len(list_id)} ASV IDs across the four consortia, and {len(list_id) - len(set_id)} ASV IDs are found more than once.")
print(f"There are a total of {len(list_seq)} ASV sequences across the four consortia, and {len(list_seq) - len(set_seq)} ASV sequences are found more than once.\n")
#print(asv_dict.items())

seq_dict = {}
df_sequence = df1["Sequence"].tolist()
df_taxa = df1["joined"].tolist()

for seq, taxa in zip(df_sequence, df_taxa):
    seq_dict[seq] = taxa
for key,value in asv_dict.items():
    asv_dict[key] = asv_dict[key], seq_dict.get(value)
    
# Dictionary is now like this: {ASV1: (seq1, taxa1), ASV2: (seq2, taxa2), ...}

# Third sanity check. Looking to see if the ASV IDs with the same sequence do point to the same taxa    
#This second sanity check confirms that there are 236 duplicates in the ASV sequences from the appended dictionary
duplicate_pairs = {}
seen_sequences = {}
for key, value in asv_dict.items():
    seq = value[0]  # Since asv_dict was updated to (seq, taxa) tuples, storing sequence
    if seq in seen_sequences:
        duplicate_pairs[seq] = (seen_sequences[seq], key)  # (original_id, duplicate_id)
    else:
        seen_sequences[seq] = key

k = 0
for (asv1,asv2) in duplicate_pairs.values():
    if(asv_dict[asv1][1] != asv_dict[asv2][1]):
        #print(f"ASV ID {asv1} and ASV ID {asv2} have the same sequence but different taxa, which is unexpected.")
        k = k + 1
    else:
        k = k
        #print(f"ASV ID {asv1} and ASV ID {asv2} have the same sequence and the same taxa, as expected.")
        
print(f"There are {k} pairs of ASV IDs with the same sequence but different taxa.")
    #print(f"ASV ID {asv1} and ASV ID {asv2} have the same sequence, as expected.")

for(key,value) in asv_dict.items():
    print(key)
    print(value)

'''

#%%
df1["Genus_species"] = (df1["genus_final"].str.strip() + " " + df1["species_final"].str.strip())
asv_dict = df1.set_index("ASV_Seq")["Genus_species"].to_dict()

for com in communities:
    for file in glob.glob(f"PICRUSt2.3/{com}_output/{com}_pathway_inference/*_strat.tsv"):
        
        df_in = pd.read_csv(file, sep = "\t")
        df_in["Taxa"] = df_in["sequence"].map(asv_dict).fillna("Unknown")
        column_to_move = df_in.pop("Taxa")
        df_in.insert(2, "Taxa", column_to_move)
        #print(df_in["sequence"].nunique())
        #print(df_in)
        #print(df_in["Taxa"])
        df_in.to_csv(f"{file.replace('_strat.tsv','_strat_with_taxa.tsv')}", sep = "\t", index = False)



# %%
