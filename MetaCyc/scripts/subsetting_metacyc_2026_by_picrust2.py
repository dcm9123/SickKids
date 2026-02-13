# Daniel Castaneda Mogollon, PhD
# February 13th, 2026
# Purpose: This script was made to filter out the MetaCyc pathways from PICRUSt2
# so I can run MinPath with the database from PICRUSt2 and compare it to the default of MinPath
# and the latest MetaCyc 2026 database

#%%
import pandas as pd
import os

#%%
master_file_2026 = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/smart_tables_parsed/Master_Metacyc_pathway_file.tsv", sep="\t")

file_name = "/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_path2rxn_struc_filt_pro.txt"
picrust2_input = open(file_name, "r")
picrust2_db_formatted = open(file_name.replace(".txt","_formatted.txt"), "w")
picrust2_pwys = []

pwy_and_reaction = {}

#icrust2_pwys = picrust2_input[0].tolist()
with picrust2_input as f:
    for line in f:
        line = line.strip()
        line = line.split("\t")
        pwy_id = line[0]
        picrust2_pwys.append(pwy_id)

for element in picrust2_pwys:
    if element in master_file_2026['Pathways'].tolist():
        pwy_and_reaction[element] = master_file_2026.loc[master_file_2026['Pathways']==element,'Reactions of pathway'].values[0]
    else:
        print(f"Error: pathway {element} is in the PICRUSt2 database but not in the MetaCyc 2026 database. Exiting now.")
        #exit()

pwy_and_reaction_df = (
    pd.DataFrame.from_dict(pwy_and_reaction, orient="index", columns=["Reactions"])
    .reset_index()
    .rename(columns={"index": "Pathway"})
)

pwy_and_reaction_df["Reactions"] = (
    pwy_and_reaction_df["Reactions"]
    .astype(str)
    .str.split(",")
)

pwy_and_reaction_df_final = (
    pwy_and_reaction_df
    .explode("Reactions")
    .assign(Reactions=lambda df: df["Reactions"].str.strip())
    .loc[lambda df: df["Reactions"].ne("")]
    .reset_index(drop=True)
)

pwy_and_reaction_df_final.to_csv(
    file_name.replace(".txt", "_pwy_rxn_exploded.tsv"),
    sep="\t",
    index=False,
)

print(pwy_and_reaction_df_final)


# %%
