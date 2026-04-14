# Daniel Castaneda Mogollon, PhD
# 03/18/2025
# Purpose: To test the difference in EC count between ASVs and the taxon it belongs to by using different penalties in the edge_exponent of PICRUSt2


#%%
from ast import Global

import pandas as pd
import os

#%%

path = "/Users/danielcm/Desktop/SickKids/hsp_testing"
os.chdir(path)

#%%
community = "s5"
type = "raw"
# Reading files and storing them as data frames
df_no_penalty = pd.read_csv(f"{type}/{community}_output/{community}_EC_nsti_predicted.tsv", sep='\t', header=0)
#df_pen_2 = pd.read_csv("ns6_EC_nsti_predicted_penalty2.tsv", sep='\t', header=0)
#df_pen_5 = pd.read_csv("ns6_EC_nsti_predicted_penalty5.tsv", sep='\t', header=0)
#df_pen_10 = pd.read_csv("ns6_EC_nsti_predicted_penalty10.tsv", sep='\t', header=0)

# Reading Anthony's file with sample ID and its corresponding genome ID
df_community_list = pd.read_csv("/Users/danielcm/Desktop/SickKids/Final_community_list_for_Daniel.csv", sep=',', header=0)

# Reading the metagenome file with EC counts for each genome ID that I made after annotating with Prokka and EggNOG-mapper
df_metagenome1 = pd.read_csv(f"{type}/EC_for_picrust2_renamed.tsv",sep='\t', header=0)

# Laura's file
#master_file = "../../../Phyloseq2/FemMicro_final_364_collapsed_all_JULY_TO_SHARE_20260306.csv"
if type == "filtered":
    master_file = f"/Users/danielcm/Desktop/SickKids/Phyloseq2/{community.upper()}_filtered_ASVs_count600_len400_prev20.csv"
else:
    master_file = f"/Users/danielcm/Desktop/SickKids/Phyloseq2/{community.upper()}_unfiltered.csv"

df_master = pd.read_csv(master_file, sep=',', header=0)


# Subsetting Laura's file

def subsetting_master_file(community):
    global df_master
    # Keeping columns of interest
    metadata_cols = ["asv_id","genus_final", "curated_species_femmicro", "asv_seq","expected_communities","asv_len"]
    # 'Grepping' by the community it belongs to, careful with S5 though... 
    sample_cols = [col for col in df_master.columns if col.startswith("plate") and (f"_{community}_") in col]
    # Sanity check
    print("There are", len(sample_cols), "columns for the", community, "community to analyze.")
    df_master_subset = df_master[metadata_cols + sample_cols].copy()
    df_master_subset["Taxa_concatenated"] = df_master_subset["genus_final"] + " " + df_master_subset["curated_species_femmicro"].copy()
    return df_master_subset

df_master_subset = subsetting_master_file(community)
print(df_master_subset["Taxa_concatenated"].head())

# The master file now has a subset of the columns of interest plus the concatenated taxa


def subsetting_rest_of_files(df_penalty):
    df_EC_calculation = df_penalty
    # Removing any unwanted tips of the tree that are not in the community of interest, careful with S5 though...
    #df_EC_calculation = df_EC_calculation[~df_EC_calculation["sequence"].str.startswith(f"S_{community}_")]

    # Subsetting Anthony's file to genome and sample ID
    df_community = df_community_list[df_community_list[["Internal sample name", "Taxonomy (GTDB-tk)"]].notna().all(axis=1)]
    df_community = df_community[["Internal sample name", "Taxonomy (GTDB-tk)"]]
    ncol = df_EC_calculation.shape[1]
    print(f"There are {len(df_EC_calculation)} ASVs and {ncol} ECs in the prediction data frame.")
    
    df_EC_calculation = df_EC_calculation.rename(columns = lambda col: col.replace("EC:", "ASV_EC:") if col.startswith("EC:") else col)
    return(df_EC_calculation, df_community)

#df_ec_calculation, df_community = subsetting_rest_of_files("NS6", df_no_penalty)
#df_ec_calculation, df_community = subsetting_rest_of_files("NS6", df_pen_2)
#df_ec_calculation, df_community = subsetting_rest_of_files("NS6", df_pen_5)
df_ec_calculation, df_community = subsetting_rest_of_files(df_no_penalty)


def merging_and_sanity_check(df_EC_calculation, df_master_subset, df_community):
    global df_metagenome1
    df_master_and_ec_prediction = df_master_subset.merge(df_EC_calculation, left_on = "asv_seq", right_on = "sequence", how = "inner")
    
    if df_master_and_ec_prediction["sequence"].equals(df_master_and_ec_prediction["asv_seq"]):
        print("The sequences match after merging.")
        df_master_and_ec_prediction["ASV_sequence_match?"] = "Yes"
    else:
        print("The sequences do not match after merging, THIS SHOULD NEVER HAPPEN.")
        df_master_and_ec_prediction["ASV_sequence_match?"] = "No"
    
        # Removing duplicates from v1v9 (which have the same profile), and renaming EC columns to Genome EC prediction 
    df_metagenome1["assembly"] = df_metagenome1["assembly"].apply(lambda x: x.replace("_v1v9",""))
    df_metagenome1["assembly"] = df_metagenome1["assembly"].apply(lambda x: "_".join(x.split("_")[0:4]))
    df_metagenome1 = df_metagenome1.rename(columns = lambda col: col.replace("EC:","Genome_EC:") if col.startswith("EC:") else col)
        
          
    df_metagenome_check = df_metagenome1.copy()
    df_metagenome_check["assembly_original"] = df_metagenome_check["assembly"]

    df_metagenome_check["assembly"] = df_metagenome_check["assembly"].apply(lambda x: x.replace("_v1v9",""))
    df_metagenome_check["assembly"] = df_metagenome_check["assembly"].apply(lambda x: "_".join(x.split("_")[0:4]))
    #print(df_metagenome_check["assembly"])
    genome_ec_cols = [col for col in df_metagenome_check.columns if col.startswith("EC:")]

    problematic_assemblies = []

    for assembly_id, group in df_metagenome_check.groupby("assembly"):
        if len(group) > 1:
            nunique_per_ec = group[genome_ec_cols].nunique()
            for item in nunique_per_ec.items():
                if item[1] > 1:
                    print(f"Assembly {assembly_id} has more than one unique value for EC column {item[0]}:")
                    print(group[["assembly_original", item[0]]])
            if (nunique_per_ec > 1).any():
                problematic_assemblies.append(assembly_id)

    print("Problematic collapsed assemblies:", len(problematic_assemblies))
    print(problematic_assemblies[:20])
    
    
    df_metagenome_f = df_metagenome1.drop_duplicates(subset=["assembly"])
    df_metagenome_f["assembly"] = df_metagenome_f["assembly"].apply(lambda x: x.replace("_v1v9",""))
    df_metagenome_f["assembly"] = df_metagenome_f["assembly"].apply(lambda x: "_".join(x.split("_")[0:4]))


    df_assembly_and_taxa = df_metagenome_f.merge(df_community, left_on = "assembly", right_on = "Internal sample name", how = "left")
    df_final = df_assembly_and_taxa.merge(df_master_and_ec_prediction, left_on = "Taxonomy (GTDB-tk)", right_on = "Taxa_concatenated", how = "inner")
    
    df_final["Taxonomy_Match?"] = df_final.apply(lambda row: "Yes" if row["Taxonomy (GTDB-tk)"] == row["Taxa_concatenated"] else "No", axis=1)
    df_final = df_final.dropna(subset = ["asv_seq"])
    asv_seq = df_final.pop("asv_seq")
    df_final = df_final.dropna(subset = ["Internal sample name"])
    internal_sample_name = df_final.pop("Internal sample name")
    asv_id = df_final.pop("asv_id")
    taxonomy_gtbtk = df_final.pop("Taxonomy (GTDB-tk)")
    genus_final = df_final.pop("genus_final")
    curated_species_femmicro = df_final.pop("curated_species_femmicro")
    expected_communities = df_final.pop("expected_communities")
    asv_len = df_final.pop("asv_len")
    taxonomy_match = df_final.pop("Taxonomy_Match?")
    metadata_NSTI = df_final.pop("metadata_NSTI")
    sequence_match = df_final.pop("ASV_sequence_match?")
    df_final.insert(1, "Internal sample name", internal_sample_name)
    df_final.insert(0, "asv_id", asv_id)
    df_final.insert(2, "Taxonomy (GTDB-tk)", taxonomy_gtbtk)
    df_final.insert(3, "genus_final", genus_final)
    df_final.insert(4, "curated_species_femmicro", curated_species_femmicro)
    df_final.insert(5, "asv_seq", asv_seq)
    df_final.insert(6, "expected_communities", expected_communities)
    df_final.insert(7, "asv_len", asv_len)
    df_final.insert(8, "Taxonomy_Match?", taxonomy_match)
    df_final.insert(9, "Metadata_NSTI", metadata_NSTI)
    df_final.insert(10, "ASV_sequence_match?", sequence_match)
    
    #df_final["EC count difference"] = df_final.apply(lambda row: row.filter(like="ASV_EC:").sum() - row.filter(like="Genome_EC:").sum(), axis=1)
    #df_final["EC count difference"] = abs(df_final["EC count difference"])
    genome_ecs = df_final.filter(like="Genome_EC:").copy()
    asv_ecs = df_final.filter(like="ASV_EC:").copy()
    
    genome_ecs.columns = genome_ecs.columns.str.replace("Genome_EC:", "", regex=False)
    asv_ecs.columns = asv_ecs.columns.str.replace("ASV_EC:", "", regex=False)

    difference = genome_ecs - asv_ecs
    difference = abs(difference)
    
    df_final["EC_total_count_difference"] = difference.sum(axis=1)
    df_final["EC_count_difference"] = (difference != 0).sum(axis=1)

    #print(difference)
    
    min_values = df_final.groupby("asv_seq")["EC_total_count_difference"].min().reset_index()
    df_final_ties = df_final.merge(min_values, on=["asv_seq", "EC_total_count_difference"])
    
    df_final_ties = df_final_ties.sort_values(["asv_seq", "Internal sample name"])
    df_final_no_ties = df_final_ties.drop_duplicates(subset=["asv_seq"], keep="first")
    
    #print(f"There were a total of {df_metagenome.shape[0]} assemblies ID, and after removing duplicates, there are {df_metagenome_f.shape[0]} unique assembly IDs in the metagenome data frame.")
    
    return(df_final_no_ties)

x =merging_and_sanity_check(df_ec_calculation, df_master_subset, df_community)  
x.to_csv(f"{type}/{type}_{community}_final_nn_no_ties.csv", index=False)    
    

# %%
