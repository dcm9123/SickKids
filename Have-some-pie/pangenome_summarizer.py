# Daniel Castaneda Mogollon, PhD
# June 16, 2026
# Purpose: This script takes the input files from the pangenome summary from multiple files from Anvi'o, and 
# creates a new table with the features of interest.

#%%
#import pip
#pip.main(['install', 'pandas'])
import pandas as pd
import os
#%%
#%%
path = "/Users/danielcm/Desktop/SickKids/Anvio/07_PANGENOMICS/"
genome1 = "AFINEGOLDII-PROJECT/"
misc_data = "SUMMARIZE/misc_data_layers/"
f_in1 = pd.read_csv(os.path.join(path,genome1, "afinegoldii_genomes.txt"), sep="\t")
f_in2 = pd.read_csv(os.path.join(path,genome1,misc_data, "default.txt"), sep="\t")
#%%

#%%
def get_genome_ids(file_in):
    num_genomes = len(file_in)
    keys = file_in["name"].tolist()
    consortium_list = file_in["contigs_db_path"].tolist()
    cons_list = []
    
    for item in consortium_list:
        item = item.split("/")[-1]
        item = item.split("_")[1]
        cons_list.append(item)
    
    genome_dict = dict(zip(keys, cons_list))
    genome_dict = dict(sorted(genome_dict.items()))
    return(num_genomes, genome_dict)
    #genome_dict["path"] = file_in["genome"].tolist()
#%%

#%%
def getting_genome_stats(file_in):
    genome_dict2 = {}
    id = file_in["layers"].tolist()
    total_length = file_in["total_length"].tolist()
    gc_content = file_in["gc_content"].tolist()
    percent_completion = file_in["percent_completion"].tolist()
    num_genes = file_in["num_genes"].tolist()
    singleton_gene_clusters = file_in["singleton_gene_clusters"].tolist()
    num_gene_clusters = file_in["num_gene_clusters"].tolist()
    genome_dict2.keys = id
    genome_dict2["total_length"] = total_length
    genome_dict2["gc_content"] = gc_content
    genome_dict2["percent_completion"] = percent_completion
    genome_dict2["num_genes"] = num_genes
    genome_dict2["singleton_gene_clusters"] = singleton_gene_clusters
    genome_dict2["num_gene_clusters"] = num_gene_clusters
    return(genome_dict2)
#%%

#%%
num_genomes,genome_ids = get_genome_ids(f_in1)

print(genome_ids)      
# %%
