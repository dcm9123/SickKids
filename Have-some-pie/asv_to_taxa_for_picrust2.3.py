# Daniel Castaneda Mogollon, PhD
# February 21st, 2026
# Purpose: To crossmap and convert each ASV from PICRUSt2.3 to its corresponding taxa from PICRUSt2.3 run

#%%
import pandas as pd
import os
import glob as glob
#%%

path = "/Users/danielcm/Desktop/SickKids/"
os.chdir(path)

asv_dict = {}

df1 = pd.read_csv("Phyloseq/t1d_db2_merged_raw_taxa_table.csv")
asv_ids = df1.iloc[:,0]

for row in df1.itertuples():
    asv_dict[row[1]] = row.genus_final,row.species_final

community = ["ns1","ns6","s2","s5"]

for com in community:
    for file in glob.glob(f"PICRUSt2/{com}_pathway_inference/*_strat.tsv"):
        df_in = pd.read_csv(file, sep = "\t")
        df_in["Taxa"] = df_in["sequence"].map(asv_dict)
        df_in["Taxa"] = df_in["Taxa"].apply(lambda x: str(x).replace("(","").replace(")","").replace(",","_").replace("'","").replace(" ",""))
        print(file_name)
        file_name = file.split("/")[-1].replace("_strat.tsv","_strat_with_taxa.tsv")
        df_in.to_csv(f"PICRUSt2/{com}_pathway_inference/{file_name}", sep="\t", index=False)
        
for com in community:
    for file in glob.glob(f"PICRUSt2/{com}_pathway_inference/*predictions.tsv"):
        df_in = pd.read_csv(file, sep = "\t")
        df_in["Taxa"] = df_in["sequence"].map(asv_dict)
        df_in["Taxa"] = df_in["Taxa"].apply(lambda x: str(x).replace("(","").replace(")","").replace(",","_").replace("'","").replace(" ",""))
        file_name = file.split("/")[-1].replace("predictions.tsv","predictions_with_taxa.tsv")
        print(file_name)
        df_in.to_csv(f"PICRUSt2/{com}_pathway_inference/{file_name}", sep="\t", index=False)
# %%
