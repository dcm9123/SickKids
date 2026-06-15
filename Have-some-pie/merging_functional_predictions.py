# Daniel Castaneda Mogollon, PhD
# 21:54
# merging_functional_predictions.py
# Purpose: This script was generated to merge the results of the EC and KO functional prediction files from
# the 4 consortia (NS1, NS6, S2, S5) into one.

#%%
import os
import pandas as pd

#%%
#This allows me to print the entire df without any issues or constraints from python
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = "/Users/danielcm/Desktop/SickKids/PICRUSt2.6/"
os.chdir(path)
consortia_list = ["ns1","s2","ns6","s5"] #Change this if needed
function_list = ["KO","EC","pathway"] #Add KO, EC or remove it
df_list_ko = []
df_list_ec = []
df_list_pathway = []
i = 0

inocula = False # Change this to False if you want to merge the sample files instead of the inoculum 12 sample files

#%%
for consortia in consortia_list:
    for function in function_list:
        print(f"Working with: {consortia}-{function}")
        if function == "KO" or function == "EC":
            if inocula == True:
                df = pd.read_csv(f"{consortia}_inocula_output/{consortia}_{function}_metagenome_out/pred_metagenome_unstrat.tsv",sep='\t',header=0) # Change this for inocula or samples
            else:
                df = pd.read_csv(f"{consortia}_output/{consortia}_{function}_metagenome_out/pred_metagenome_unstrat.tsv",sep='\t',header=0) # Change this for inocula or samples
            df = df.set_index('function')
            if function == "KO":
                df_list_ko.append(df)
                i = i+1
            elif function == "EC":
                df_list_ec.append(df)
                i = i+1
            else:
                print(f"This should never happen. Invalid function type: {function}.")
        else:
            if inocula == True:
                df =pd.read_csv(f"{consortia}_inocula_output/{consortia}_pathway_inference/path_abun_unstrat.tsv",sep='\t',header=0) # Change this for inocula or samples
            else:
                df =pd.read_csv(f"{consortia}_output/{consortia}_pathway_inference/path_only_abun_strat_no_taxa_aggregated.tsv",sep='\t',header=0) # Change this for inocula or samples or file name
            #df = df.set_index('pathway') #Change to pathway_taxa or pathway when working with stratified pathway output vs not
            #df["pathway_taxa"] = df["pathway"] + "_" + df["Taxa"] #Change this if working with stratified pathway output vs not
            #df = df.set_index('pathway_taxa') #Change this if working with stratified pathway output vs not
            #print(len(df.index))
            #print(len(df.index.unique()))
            df = df.set_index('pathway') #Change this if working with stratified pathway output vs not
            df_list_pathway.append(df)
            i = i+1

df_ko_merged = pd.concat(df_list_ko, axis=1).fillna(0)
df_ec_merged = pd.concat(df_list_ec, axis=1).fillna(0)
df_pathway_merged = pd.concat(df_list_pathway, axis=1).fillna(0)

df_ko_merged_t = df_ko_merged.T
df_ec_merged_t = df_ec_merged.T
df_pathway_merged_t = df_pathway_merged.T

if inocula == True:
    df_ko_merged_t.to_csv("KO_inocula_merged_metagenome.tsv",index=True, sep="\t")
    df_ec_merged_t.to_csv("EC_inocula_merged_metagenome.tsv",index=True, sep="\t")
    df_pathway_merged_t.to_csv("Pathway_inocula_merged_metagenome.tsv",index=True, sep="\t")
else:
    #df_ko_merged.to_csv("KO_merged_metagenome.tsv",index=True, sep="\t")
    #df_ec_merged.to_csv("EC_merged_metagenome.tsv",index=True, sep="\t")
    #df_pathway_merged.to_csv("Pathway_merged_metagenome_strat.tsv",index=True, sep="\t")
    df_pathway_merged.to_csv("Pathway_merged_metagenome_strat.tsv",index=True, sep="\t")


for list_in in df_list_ec,df_list_pathway: # add df_list_ko if needed
    for dframe in list_in:
        print(f"The dimension of the df are: {dframe.shape}")

print(f"The dimensions of the merged KO df are {df_ko_merged_t.shape}")
print(f"The dimensions of the merged EC df are {df_ec_merged_t.shape}")
print(f"The dimensions of the merged Pathway df are {df_pathway_merged_t.shape}")


for item in df_ko_merged.index: # Check that all items in the merged df are present in at least one of the input dfs
    found = any(item in df.index for df in df_list_ko) # Check if item is in any of the input dataframes
    if not found:
        print(f"{item} not found in any KO input dataframe (unexpected)") # This should never happen

for dframe in df_list_ko:
    for value in dframe.index:
        if value not in df_ko_merged.index:
            print(f"{value} not found in the merged KO df (unexpected)") # This should never happen

for item in df_ec_merged.index:
    found = any(item in df.index for df in df_list_ec)
    if not found:
        print(f"{item} not found in any EC input dataframe (unexpected)") # This should never happen

for dframe in df_list_ec:
    for value in dframe.index:
        if value not in df_ec_merged.index:
            print(f"{value} not found in the merged EC df (unexpected)") # This should never happen

for item in df_pathway_merged.index:
    found = any(item in df.index for df in df_list_pathway)
    if not found:
        print(f"{item} not found in any pathway input dataframe (unexpected)") # This should never happen

for dframe in df_list_pathway:
    for value in dframe.index:
        if value not in df_pathway_merged.index:
            print(f"{value} not found in the merged pathway df (unexpected)") # This should never happen


# %%
