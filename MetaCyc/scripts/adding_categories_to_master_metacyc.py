# Daniel Castaneda Mogollon, PhD
# February 22nd, 2026
# Purpose: To add the category of each pathway to the Master MetaCyc database created previously. 
# using the existing metacyc curated categories file.

#%%
import pandas as pd
import os

path = "/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/"
os.chdir(path)

master_df = pd.read_csv("Master_Metacyc_pathway_file.tsv", sep="\t")
df_categories = pd.read_csv("metacyc_category_curation_2026.txt", sep="\t")

print(master_df.columns)
print(df_categories.columns)

df_categories["Unique"].is_unique #True, sanity check done

df_final = master_df.merge(df_categories.drop_duplicates("Unique"), left_on="Ontology - pathway type", right_on="Unique", how="left",
                           validate="many_to_one")

df_final.to_csv("Master_Metacyc_pathway_file_with_categories.tsv", sep="\t", index=False)

#master_df["Level1"] = master_df["Ontology - pathway type"].map(df_categories)

# %%
