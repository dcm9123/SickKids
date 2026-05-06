# Improvement from Have-some-pie


# Loading packages ----------------------------------------------------------
loading_packages = function() {
    library(igraph)
    library(compositions)
    library(RColorBrewer)
    library(gplots)
    library(beepr)
    library(ggraph)
    library(psych)
    library(tidyr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(tidygraph)
    library(ggrepel)
    library(RCy3)
    print("All packages loaded successfully")
}

loading_packages()

working_path = "/Users/danielcm/Desktop/SickKids/"
setwd(working_path)


pathway_file = "PICRUSt2.6/pathway_merged_metagenome.tsv"
taxa_file = "Phyloseq2/NS1_filtered_ASVs_count600_len400_prev20.csv"
maaslin_file = "Maaslin2.6/w9w10/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_Consortium,S2_w9w10/NS1_w9w10_vs_S2_w9w10_enriched_pathways_with_categories.csv"
metadata_file = "Metadata/Danska_diabetes_metadata364_20260409.csv"
output_path = "T1D/Maaslin2.6/Networks/NS1"
filtering_category = "Week_and_Consortium"
filtering_value = "NS1_w9w10"
alpha_sign = 0.05
rho_minimum = 0.7

tmp = read.csv(taxa_file, header = TRUE, row.names = 1)
pathw = read.csv(pathway_file, header = TRUE, row.names = 1, sep = "\t")
mtd = read.csv(metadata_file, header = TRUE, row.names = 1)

rownames(pathw)
colnames(tmp)
rownames(mtd)

# INPUT:
# This function takes as input a taxa file for a consortium (NS1, NS6, etc) with genus_final, curated_species_femmicro, 
# and the columns having sample names, and creates a new column with the combined genus and species name, 
# and then collapses the table to genus_species level by summing the counts for ASVs that belong to the same genus_species. 
# It also reads in the pathway file and makes sure the sample names match between the taxa and pathway tables, and if not, identifies which ones are missing from each table. It returns a list with the processed pathway dataframe and the processed taxa dataframe.

# OUTPUT: Returns a count table with the names matching the pathway_merged_metagenome file from PICRUSt2
reading_files_and_matching_names = function(pathway_file, taxa_file, metadata_file){
    # Reading pathway, taxa, and metadata files
    pathway_df = as.data.frame(read.table(pathway_file, header = TRUE, sep = "\t", row.names = 1))
    taxa_df = as.data.frame(read.csv(taxa_file, header = TRUE, sep = ",", row.names = 1))
    metadata_df = read.csv(metadata_file, header = TRUE, sep = ",", row.names = 1) # Read with row.names=1 so Id values become rownames

    # Making a new genus and species combined from the matrix
    taxa_df$genus_species = paste(taxa_df$genus_final, taxa_df$curated_species_femmicro, sep = "_")

    # Identifying columns that are not sample columns, and in another variable, keep only sample columns.
    non_asv_columns = c("genus_final","curated_species_femmicro","expected_communities","asv_len","asv_seq","sum","nonzero_count","genus_species")
    sample_columns = setdiff(colnames(taxa_df), non_asv_columns)

    # Collapsing the taxa table to genus_species level by summing the counts for ASVs that belong to the same genus_species
    taxa_df_collapsed = taxa_df %>%
        group_by(genus_species) %>%
        summarise(across(all_of(sample_columns), sum, na.rm = TRUE))
    # Store genus_species names before any column filtering (tibbles don't support rownames, assigned at the end)
    genus_species_names = gsub(" ", "_", taxa_df_collapsed$genus_species)
    rownames(taxa_df_collapsed) = genus_species_names
    # Getting the right nomenclature
    colnames(pathway_df) = gsub("^X", "", colnames(pathway_df)) # Removing any leading X from column names if present
    colnames(pathway_df) = gsub("\\.", "-", colnames(pathway_df)) # Replacing dots with dashes to match metacyc names

    print(paste0("Originally, there are a total of ", ncol(taxa_df_collapsed) - 1, " samples in the taxa table and ", ncol(pathway_df), " samples in the pathway table.")) # -1 to exclude the genus_species column
 
    # drop the genus_species column since we will assign it as rownames at the end, and we don't want it to interfere with the matching of sample names between the taxa and pathway tables
   
    # Replace the names of the samples in the taxa file with the updated SRA names
    # setNames(new_names, old_names): rownames(metadata_df) are the Id values (current taxa column names), metadata_df$SRA_sample_name are the target names

    mapping <- setNames(metadata_df$SRA_sample_name, rownames(metadata_df))

    taxa_df_collapsed_renamed <- taxa_df_collapsed %>%
        rename(any_of(mapping))

    pattern = "(?i)(inoc|inocula|inoculum|control|ctrl|positive|negative|defined_inocula)"    # Keeping only real sample columns, dropping controls/inocula
    taxa_df_collapsed_f = taxa_df_collapsed_renamed %>%
        select(-matches(pattern, perl = TRUE))
    taxa_df_collapsed_f = taxa_df_collapsed_f %>%
        select(-genus_species) # drop the genus_species column since we will assign it as rownames at the end, and we don't want it to interfere with the matching of sample names between the taxa and pathway tables

    dropped_cols = setdiff(colnames(taxa_df_collapsed_renamed), c("genus_species", colnames(taxa_df_collapsed_f)))
    print(paste0("After removing controls and inocula, there are a total of ", ncol(taxa_df_collapsed_f), " samples in the taxa table."))
    print(paste0("The dropped non-plate columns are: ", paste(dropped_cols, collapse = ", ")))
 
    # Sanity check: how many taxa columns now match pathway columns
    matched = intersect(colnames(taxa_df_collapsed_f), rownames(pathway_df)) # Rownames for pathway_df are samples, and colnames for taxa
    missing_from_pathway = setdiff(colnames(taxa_df_collapsed_f), rownames(pathway_df))
    missing_from_taxa = setdiff(rownames(pathway_df), colnames(taxa_df_collapsed_f))
    print(paste0(length(matched), "/", ncol(taxa_df_collapsed_f), " taxa sample columns match pathway table rows after renaming."))
    
    if(length(missing_from_pathway) > 0){
        print(paste0("In taxa but NOT in pathway table: ", paste(missing_from_pathway, collapse = ", ")))
    }
    if(length(missing_from_taxa) > 0){
        print(paste0("In pathway table but NOT in taxa: ", paste(missing_from_taxa, collapse = ", ")))
    }

    # Assign rownames (genus_species) now that all transformations are done; must convert tibble to data.frame first
    taxa_df_final = as.data.frame(taxa_df_collapsed_f)
    print(colnames(taxa_df_final))
    rownames(taxa_df_final) = genus_species_names

    return(list(pathway_df = pathway_df, taxa_df = taxa_df_final))
}


result = reading_files_and_matching_names(pathway_file, taxa_file, metadata_file)
pathway_file2 = result$pathway_df  # Rows are pathways, columns are samples
taxa_file2 = result$taxa_df     # Rows are genus_species, columns are samples

maaslin_df = read.csv(maaslin_file, header = TRUE, sep = ",")
rownames(maaslin_df) = maaslin_df$feature

# Sanity check: Make sure all the pathways from maaslin2 are present in the merged_pathway file from PICRUSt2
if(all(rownames(maaslin_df) %in% colnames(pathway_file2))){
    print("All features in the maaslin file are present in the pathway file.")
} else {
    missing_in_pathway = setdiff(rownames(maaslin_df), colnames(pathway_file2))
    stop(paste0("The following features in the maaslin file are NOT present in the pathway file: ", paste(missing_in_pathway, collapse = ", ")))
}

pathway_file2_f = pathway_file2 %>%
    select(all_of(rownames(maaslin_df))) # keep only the pathways that are in the maaslin file, and in the same order, so we can easily match the maaslin results to the pathway abundances when we get to the network construction step

pathway_file2_f_t = t(pathway_file2_f) # transpose so rows are samples and columns are pathways, to match the taxa file format and make it easier to combine them later for the network construction step
pathway_file2_f_t = as.data.frame(pathway_file2_f_t)


pathway_file2_f_t = pathway_file2_f_t %>%
    select(all_of(colnames(taxa_file2))) # keep only the samples that are in the taxa file, and in the same order, so we can easily match the taxa abundances to the pathway abundances when we get to the network construction step

colnames(taxa_file2) # should be the same as rownames(pathway_file2_f_t)
View(taxa_file2)

ncol(pathway_file2_f) # should be the same as nrow(maaslin_df)






normalizing_and_filtering = function(pathway_file, taxa_file, maaslin_filt_file, metadata_file, output_path, output_file_name_net, network_categories_output, 
filtering_category, filtering_value, alpha_sign, rho_minimum, ref_value){
    
    # Reading master files as tables and matrices
    maaslin_filt = read.table(maaslin_filt_file, header = TRUE, sep = "\t")
    metadata_df = read.csv(metadata_file, header = TRUE, sep = ",", row.names = 1)
    pathway_mat = as.matrix(read.table(pathway_file, header = TRUE, sep = "\t", row.names = 1))
    taxa_mat = as.matrix(read.csv(taxa_file, header = TRUE, sep = ",", row.names = 1))
    
    # Renaming columns, keeping it consistent
    colnames(pathway_mat) = gsub("^X", "", colnames(pathway_mat)) # Removing any leading X from column names if present
    colnames(pathway_mat) = gsub("\\.", "-", colnames(pathway_mat)) # Replacing dots with dashes to match metacyc names




}

