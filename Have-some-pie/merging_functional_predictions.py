# Daniel Castaneda Mogollon, PhD
# 21:54
# merging_functional_predictions.py
# Purpose: This script was generated to merge the results of the EC and KO functional prediction files from
# the 4 consortia (NS1, NS6, S2, S5) into one.

import pandas as pd
from pathlib import Path

#%%
#This allows me to print the entire df without any issues or constraints from python
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = Path("/Users/danielcm/Desktop/SickKids/PICRUSt2.6/")
consortia_list = ["ns1","s2","ns6","s5"] #Change this if needed
function_list = ["KO","EC","pathway"] #Add KO, EC or remove it
pathway_file_name = "path_abun_unstrat.tsv" # Change to path_abun_strat_with_taxa_aggregated.tsv, etc. if needed
valid_functions = {"KO", "EC", "pathway"}
df_list_ko = []
df_list_ec = []
df_list_pathway = []

inocula = False # Change this to False if you want to merge the sample files instead of the inoculum 12 sample files


def get_input_file(consortia, function, inocula):
    if function in ["KO", "EC"]:
        output_type = "inocula_output" if inocula else "output"
        return path / f"{consortia}_{output_type}" / f"{consortia}_{function}_metagenome_out" / "pred_metagenome_unstrat.tsv"

    output_type = "inocula_output" if inocula else "output"
    return path / f"{consortia}_{output_type}" / f"{consortia}_pathway_inference" / pathway_file_name


def prepare_dataframe(file_to_read, function):
    df = pd.read_csv(file_to_read, sep="\t", header=0)

    if function in ["KO", "EC"]:
        index_cols = ["function"]
    elif {"pathway", "sequence", "Taxa"}.issubset(df.columns):
        index_cols = ["pathway", "sequence", "Taxa"]
    elif "pathway_taxa" in df.columns:
        index_cols = ["pathway_taxa"]
    elif {"pathway", "Taxa"}.issubset(df.columns):
        index_cols = ["pathway", "Taxa"]
    elif {"pathway", "sequence"}.issubset(df.columns):
        index_cols = ["pathway", "sequence"]
    else:
        index_cols = ["pathway"]

    missing_cols = [col for col in index_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"{file_to_read} is missing expected ID column(s): {missing_cols}")

    sample_cols = [col for col in df.columns if col not in index_cols]
    duplicated_samples = df[sample_cols].columns[df[sample_cols].columns.duplicated()].tolist()
    if duplicated_samples:
        raise ValueError(f"{file_to_read} has duplicated sample column(s): {duplicated_samples}")

    df = df.set_index(index_cols)
    if df.index.has_duplicates:
        duplicated_ids = df.index[df.index.duplicated()].unique().tolist()[:10]
        raise ValueError(f"{file_to_read} has duplicated feature ID(s), first examples: {duplicated_ids}")

    return df


def merge_dataframes(df_list, function):
    if not df_list:
        raise ValueError(f"No {function} dataframes were loaded. Check function_list and input paths.")
    merged = pd.concat(df_list, axis=1).fillna(0)
    duplicated_samples = merged.columns[merged.columns.duplicated()].tolist()
    if duplicated_samples:
        raise ValueError(f"Merged {function} dataframe has duplicated sample column(s): {duplicated_samples}")
    return merged


#%%
invalid_functions = set(function_list) - valid_functions
if invalid_functions:
    raise ValueError(f"Invalid function type(s): {sorted(invalid_functions)}")

for consortia in consortia_list:
    for function in function_list:
        print(f"Working with: {consortia}-{function}")
        file_to_read = get_input_file(consortia, function, inocula)
        df = prepare_dataframe(file_to_read, function)

        if function == "KO":
            df_list_ko.append(df)
        elif function == "EC":
            df_list_ec.append(df)
        elif function == "pathway":
            df_list_pathway.append(df)
        else:
            raise ValueError(f"Invalid function type: {function}")

df_ko_merged = merge_dataframes(df_list_ko, "KO") if "KO" in function_list else None
df_ec_merged = merge_dataframes(df_list_ec, "EC") if "EC" in function_list else None
df_pathway_merged = merge_dataframes(df_list_pathway, "pathway") if "pathway" in function_list else None

df_ko_merged_t = df_ko_merged.T if df_ko_merged is not None else None
df_ec_merged_t = df_ec_merged.T if df_ec_merged is not None else None
df_pathway_merged_t = df_pathway_merged.T if df_pathway_merged is not None else None

if inocula == True:
    if df_ko_merged_t is not None:
        df_ko_merged_t.to_csv(path / "KO_inocula_merged_metagenome.tsv",index=True, sep="\t")
    if df_ec_merged_t is not None:
        df_ec_merged_t.to_csv(path / "EC_inocula_merged_metagenome.tsv",index=True, sep="\t")
    if df_pathway_merged_t is not None:
        df_pathway_merged_t.to_csv(path / "Pathway_inocula_merged_metagenome.tsv",index=True, sep="\t")
else:
    if df_ko_merged is not None:
        df_ko_merged.to_csv(path / "KO_merged_metagenome.tsv",index=True, sep="\t")
    if df_ec_merged is not None:
        df_ec_merged.to_csv(path / "EC_merged_metagenome.tsv",index=True, sep="\t")
    if df_pathway_merged is not None:
        #df_pathway_merged.to_csv("Pathway_merged_metagenome_strat.tsv",index=True, sep="\t")
        df_pathway_merged.to_csv(path / "Pathway_merged_metagenome.tsv",index=True, sep="\t")


for list_in in df_list_ec,df_list_pathway: # add df_list_ko if needed
    for dframe in list_in:
        print(f"The dimension of the df are: {dframe.shape}")

if df_ko_merged_t is not None:
    print(f"The dimensions of the merged KO df are {df_ko_merged_t.shape}")
if df_ec_merged_t is not None:
    print(f"The dimensions of the merged EC df are {df_ec_merged_t.shape}")
if df_pathway_merged_t is not None:
    print(f"The dimensions of the merged Pathway df are {df_pathway_merged_t.shape}")


if df_ko_merged is not None:
    for item in df_ko_merged.index: # Check that all items in the merged df are present in at least one of the input dfs
        found = any(item in df.index for df in df_list_ko) # Check if item is in any of the input dataframes
        if not found:
            print(f"{item} not found in any KO input dataframe (unexpected)") # This should never happen

    for dframe in df_list_ko:
        for value in dframe.index:
            if value not in df_ko_merged.index:
                print(f"{value} not found in the merged KO df (unexpected)") # This should never happen

if df_ec_merged is not None:
    for item in df_ec_merged.index:
        found = any(item in df.index for df in df_list_ec)
        if not found:
            print(f"{item} not found in any EC input dataframe (unexpected)") # This should never happen

    for dframe in df_list_ec:
        for value in dframe.index:
            if value not in df_ec_merged.index:
                print(f"{value} not found in the merged EC df (unexpected)") # This should never happen

if df_pathway_merged is not None:
    for item in df_pathway_merged.index:
        found = any(item in df.index for df in df_list_pathway)
        if not found:
            print(f"{item} not found in any pathway input dataframe (unexpected)") # This should never happen

    for dframe in df_list_pathway:
        for value in dframe.index:
            if value not in df_pathway_merged.index:
                print(f"{value} not found in the merged pathway df (unexpected)") # This should never happen


# %%
