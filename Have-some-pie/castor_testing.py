# Daniel Castaneda Mogollon, PhD
# 03/18/2025
# Purpose: To test the difference in EC count between ASVs and the taxon it belongs to by using different penalties in the edge_exponent of PICRUSt2


#%%
import pandas as pd
import os

#%%
path = "/Users/danielcm/Desktop/SickKids/PICRUSt2.3/ns6_output/castor_testing/"
os.chdir(path)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
#%%
# Reading files and storing them as data frames
df_no_penalty = pd.read_csv("ns6_EC_nsti.predicted.tsv", sep='\t', header=0)
df_pen_2 = pd.read_csv("ns6_EC_nsti_predicted_penalty2.tsv", sep='\t', header=0)
df_pen_5 = pd.read_csv("ns6_EC_nsti_predicted_penalty5.tsv", sep='\t', header=0)
df_pen_10 = pd.read_csv("ns6_EC_nsti_predicted_penalty10.tsv", sep='\t', header=0)
df_community_list = pd.read_csv("../../../Final_community_list_for_Daniel.csv", sep=',', header=0)

df_metagenome = pd.read_csv("../../trait_tables/EC_for_picrust2_renamed.tsv", sep='\t', header=0)

master_file = "../../../Phyloseq2/FemMicro_final_364_collapsed_all_JULY_TO_SHARE_20260306.csv"
df_master = pd.read_csv(master_file, sep=',', header=0)

#%%
df_master = df_master[["asv_id","genus_final","curated_species_femmicro","asv_seq"]]
df_no_penalty = df_no_penalty[~df_no_penalty["sequence"].str.startswith("S_NS6_")]
df_no_penalty_analysis = df_no_penalty.merge(df_master, left_on = "sequence", right_on = "asv_seq", how = "left")

if df_no_penalty_analysis["sequence"].equals(df_no_penalty_analysis["asv_seq"]):
    df_no_penalty_analysis["Match"] = "Yes"
    
df_no_penalty_analysis["Final_taxa"] = df_no_penalty_analysis["genus_final"] + " " + df_no_penalty_analysis["curated_species_femmicro"]

df_community_list = df_community_list[["Internal sample name","Taxonomy (GTDB-tk)"]]

df_metagenome_with_taxa = df_metagenome.merge(df_community_list, left_on = "assembly", right_on = "Internal sample name", how = "left")
# %%
meta_cols_metagenome = ["assembly", "Internal sample name", "Taxonomy (GTDB-tk)"]
meta_cols_asv = ["sequence", "metadata_NSTI", "asv_id", "genus_final", 
                 "curated_species_femmicro", "asv_seq", "Match", "Final_taxa"]


ec_columns = [col for col in df_metagenome_with_taxa.columns if col not in meta_cols_metagenome]
print("There are", len(ec_columns), "ECs in the metagenome data frame.")

genome_ec = df_metagenome_with_taxa.groupby("Taxonomy (GTDB-tk)")[ec_columns].sum()
print("There are", len(genome_ec), "unique taxa in the metagenome data frame.")

asv_ec = df_no_penalty_analysis.copy()
asv_ec = asv_ec[~asv_ec["Final_taxa"].isna()]
print("There are", len(asv_ec), "ASVs with assigned taxa in the ASV data frame.")

merged = asv_ec.merge(
    genome_ec.reset_index(),
    left_on="Final_taxa",
    right_on="Taxonomy (GTDB-tk)",
    how="left",
    suffixes=("_asv", "_genome")
)

print("There are", len(merged), "rows in the merged data frame.")

# %%
