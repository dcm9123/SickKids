# Daniel Castaneda Mogollon, PhD
# October 13th, 2025
# Purpose: To create igraph networks from microbiome data

# install.packages("igraph")
# install.packages("ggraph")
# install.packages("psych")
# BiocManager::install("RCy3")


welcome_function = function() {
    print("Welcome to the Diammatics T1D Network Analysis Script!")
    print("This script will help you create and visualize networks based on microbiome taxa and pathway data.")
    print("Make sure to have all required packages installed and your data files ready, which can be found in the loading packages function")
    print("---------------------------------------------------")
    print("INPUT FILES AND PATHS NEEDED: ")
    print("1. Pathway file: A pathway file having all the pathways for all your groups, timepoints, etc. In my case it is the pathways of all communities")
    print(" NS1, S2, S5, NS6 merged together from PICRUSt2 predictions. It should be a .tsv file.")
    print("2. Taxa file: A taxa file for a particular group, i.e. ps_ns1_final.csv for NS1 community. It should be a .csv file.")
    print("3. Metadata file: A metadata file containing all the sample information, including the filtering category and value. It should be a .csv file.")
    print("4. Maasline filtered file: A maaslin2 output file that has the results from a pairwise comparison between two groups, i.e. NS1_w5w6 vs NS1_w9w10. It should be a .tsv file.")
    print("5. Output path: A path to a folder where you want to save the output files, i.e. /Users/danielcm/Desktop/diammatics/T1D/Networks/NS1_w5w6_vs_NS1_w9w10")
    print("INPUT VALUES NEEDED: ")
    print("1. Filtering category: A category that you wish to plot and analyze, as part of a column in your metadata file, i.e. Week_and_consortia")
    print("2. Filtering value: The value within the filtering category that you wish to analyze, i.e. NS1_w9w10")
    print("3. An alpha significance value: This is for plotting significant correlations, i.e. 0.05, 0.001")
    print("4. A rho minimum value: This is for plotting edges in the network, i.e. 0.7, 0.8 (We recommend not less than 0.7 for busy graphs)")
    print("5. Community name: The name of the community you are individually analyzing in a single plot, i.e. NS1, S2, S5, NS6") 
    print("6. Timepoint: The timepoint you are analyzing, i.e. w5w6, w9w10")
    print("---------------------------------------------------")
    print("OUTPUT FILES AND FIGURES: ")
    print("1. summary_analysis.txt: A summary text file containing:
    - Number of pathways from the original-merged pathway file
    - Number of input pathways ran by Maaslin2
    - Number of significant pathways found by Maaslin2 after q-val of 0.05 and Log2FC of > |1|
    - Number of taxa of the original community file
    - Sanity check of random samples between taxa and pathway matrices, ensuring the order is the same
    - Number and name of taxa filtered by not having any mathematical value by Spearman Rho after CLR normalization
    - Number and name of pathways filtered by not having any mathematical value by Spearman Rho after CLR normalization
    - Total number of correlations calculated between taxa and pathways
    - Number of significant correlations found after applying the alpha significance value and rho minimum value
    - The range of the Spearman Rho values found in the significant correlations")
    print("2. spearman_clear_correlation_test.tsv: A table containing all the Spearman correlation values and p-values between all CLR-normalized (separately) pathways and taxa and NAs removed")
    print("3. spearman_p_val_matrix.tsv: A matrix file containing all the p-values from the Spearman correlation test between all CLR-normalized (separately) pathways and taxa and NAs removed")
    print("4. spearman_significant_matrix.tsv: A matrix file containing only the rho values of significant correlations based on the alpha significance value and rho minimum value provided")
    print("5. spearman_maaslin2_enriched_clean_correlation_test.tsv: A table containing only the Spearman correlation values and p-values between CLR-normalized (separately) pathways and taxa that were significantly enriched by maaslin2 and NAs removed")
    print("5. maaslin2_significant_correlation_heatmap_plot.png: The heatmap plot for the maaslin2 significantly enriched pathways in the community of choice")
    print("6. network_summary_analysis.txt: A summary text file containing all the relevant information about the network analysis")
    print("7. Cytoscape network session: A Cytoscape session file containing the network created, which can be further modified and analyzed in Cytoscape")
    print("---------------------------------------------------")
}

fetching_files = function(){
    
    #pathway_file = readline(prompt = "Enter the full path to the pathway file e.g., /path/to/pathway_file.tsv: ")
    #if(!file.exists(pathway_file)){
    #    stop("Pathway file does not exist. Please check the path and try again.")
    #}
    #taxa_file = readline(prompt = "Enter the full path to the taxa file e.g., /path/to/taxa_file.csv: ")
    #if(!file.exists(taxa_file)){
    #    stop("Taxa file does not exist. Please check the path and try again.")
    #}
    #maaslin_filt_file = readline(prompt = "Enter the full path to the maaslin filtered file e.g., /path/to/maaslin_filt/all_results.tsv: ")
    #if(!file.exists(maaslin_filt_file)){
    #    stop("Maaslin filtered file does not exist. Please check the path and try again.")
    #}
    #metadata_file = readline(prompt = "Enter the full path to the metadata file e.g., /path/to/metadata_file.csv: ")
    #if(!file.exists(metadata_file)){
    #    stop("Metadata file does not exist. Please check the path and try again.")
    #}
    #output_path = readline(prompt = "Enter the full path to the output directory e.g., /path/to/output_directory/: ")
    #if(!dir.exists(output_path)){
    #    print("The output directory does not exist. Creating it now...")
    #    dir.create(output_path, recursive = TRUE)
    #}

    #filtering_category = readline(prompt = "Enter the filtering category from the metadata file e.g., Week_and_consortia: ")
    #filtering_value = readline(prompt = "Enter the filtering value from the filtering category e.g., NS1_w5w6: ")
    #alpha_sign = as.numeric(readline(prompt = "Enter the alpha significance value (to be used after BH adjustment) e.g., 0.05: "))
    #rho_minimum = as.numeric(readline(prompt = "Enter the rho minimum value e.g., 0.7: "))
    #community = readline(prompt = "Enter the community name e.g., NS1: ")
    #timepoint = readline(prompt = "Enter the timepoint e.g., w5w6: ")

    
    pathway_file = "/Users/danielcm/Desktop/SickKids6PICRUSt2.6/pathway_merged_metagenome.tsv"
    taxa_file = "/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ps_ns1_final.csv"
    maaslin_filt_file = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/same_community/Pathway/maaslin2_Week_and_consortia_pathway_NS1_w5_vs_NS1_w9w10_ref_Week_and_consortia,NS1_w9w10/all_results.tsv"
    metadata_file = "/Users/danielcm/Desktop/diammatics/T1D/metadata_ps.csv"
    output_path = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/NS1"
    filtering_category = "Week_and_consortia"
    filtering_value = "NS1_w9w10"
    alpha_sign = 0.05
    rho_minimum = 0.7

    #clean_correlation_matrix_file = as.matrix(read.table("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/s2_network_output/spearman_maaslin2_enriched_clean_correlation_test.tsv", header = TRUE, sep = "\t"))

    ref_value = sub(".*,(.*?)\\/.*", "\\1", maaslin_filt_file) # Extracting the reference value from the maaslin filtered file path)

    return(c(pathway_file = pathway_file, taxa_file = taxa_file, maaslin_filt_file = maaslin_filt_file, metadata_file = metadata_file, 
                output_path = output_path, filtering_category = filtering_category, filtering_value = filtering_value, alpha_sign = alpha_sign, 
                rho_minimum = rho_minimum, community = strsplit(filtering_value,"_",)[[1]][1], timepoint = strsplit(filtering_value,"_",)[[1]][2], ref_value = ref_value))
}

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

correcting_names = function(maaslin_pathways, significant, filtering_category, filtering_value, ref_value) {
    if (significant == "yes") {
        if(filtering_value == ref_value){
            maaslin_pathways = maaslin_pathways[maaslin_pathways$qval < 0.05 & maaslin_pathways$coef < -1, ] # If the filtering value is the reference, we look for negative coefficients
            if(nrow(maaslin_pathways) == 0){
                return(NULL)
            }
        } else {
            maaslin_pathways = maaslin_pathways[maaslin_pathways$qval < 0.05 & maaslin_pathways$coef > 1, ] # If the filtering value is not the reference, we look for positive coefficients
            if(nrow(maaslin_pathways) == 0){
                return(NULL)
            }
        }
    }
    corrected_pathways_names = gsub("\\.", "-", maaslin_pathways$feature) # Replacing dots with dashes to match metacyc names
    if(length(corrected_pathways_names) == 0){
        return(NULL)
    }
    for (i in 1:nrow(maaslin_pathways)) {
        if (grepl("^X", corrected_pathways_names[i])) { # If it starts with an X, remove it
            corrected_pathways_names[i] = substr(corrected_pathways_names[i], 2, nchar(corrected_pathways_names[i])) # Remove the X at the start by taking the second character to the last one.
        }
    }
    beep(1)
    return(corrected_pathways_names)
}

sample_size_revision = function(common_samples, taxa_mat, pathway_mat){
    if (length(common_samples) < 3) {
    stop("Not enough common samples (", length(common_samples),
            ") between taxa and pathway matrices for correlation. Exiting.")
    }
    else{
        print(paste0("Number of common samples between taxa and pathway matrices for correlation: ", length(common_samples)))
    }

    if (nrow(taxa_mat) < 2 || ncol(taxa_mat) < 2) {
        message("Not enough taxa or samples for CLR in taxa_mat for ", filtering_value,
            " (dim = ", paste(dim(taxa_mat), collapse = " x "), "). Skipping.")
        return(NULL)
    }

    if (nrow(pathway_mat) < 2 || ncol(pathway_mat) < 2) {
        message("Not enough pathways or samples for CLR in pathway_mat for ", filtering_value,
            " (dim = ", paste(dim(pathway_mat), collapse = " x "), "). Skipping.")
        return(NULL)
    }
    print("Sufficient samples and features for CLR normalization and correlation analysis. Proceeding...")
    return(TRUE)
}

self_normalizing_and_filtering = function(pathway_file, taxa_file, metadata_file, output_path, filtering_category, filtering_value, alpha_sign, rho_minimum, which_correlation){
    # Reading files...
    metadata_df = read.csv(metadata_file, header = TRUE, sep = ",", row.names = 1)
    pathway_mat = as.matrix(read.table(pathway_file, header = TRUE, sep = "\t", row.names = 1))
    taxa_mat = as.matrix(read.csv(taxa_file, header = TRUE, sep = ",", row.names = 1))
    # Renaming columns...
    colnames(pathway_mat) = gsub("^X", "", colnames(pathway_mat)) # Removing any leading X from column names if present
    colnames(pathway_mat) = gsub("\\.", "-", colnames(pathway_mat)) # Replacing dots with dashes to match metacyc names
    rownames(taxa_mat) = gsub(" ", "_", rownames(taxa_mat)) # Replacing spaces with underscores in taxa names
    pathway_mat = t(pathway_mat)

    filtered_IDs = rownames(metadata_df[metadata_df[[filtering_category]] == filtering_value, , drop = FALSE])
    pathway_mat  = pathway_mat[, colnames(pathway_mat) %in% filtered_IDs, drop = FALSE]
    taxa_mat = taxa_mat[,colnames(taxa_mat) %in% filtered_IDs, drop = FALSE]

    pathway_mat = pathway_mat[rowSums(pathway_mat) > 0, ] # Remove pathways with all zeros
    taxa_mat = taxa_mat[rowSums(taxa_mat) > 0, ] # Remove taxa with all zeros

    common_samples = intersect(colnames(pathway_mat), colnames(taxa_mat)) # This wll print the common samples between the two matrices, which in this case should be 100% identical
    unique_samples = setdiff(colnames(pathway_mat), colnames(taxa_mat)) # This will print the samples that are unique to clr_pathway, in theory only week 7 should be different
    
    sample_sanity_check = sample_size_revision(common_samples, taxa_mat, pathway_mat)
    if(sample_sanity_check == TRUE){
        print("Proceeding with CLR normalization and correlation analysis...")
    }
    # CLR normalization
    clr_pathway = t(as.matrix(clr(t(pathway_mat) + 1))) # Adding a pseudocount of 1 to avoid issues with zeros
    clr_taxa = t(as.matrix(clr(t(taxa_mat) + 1))) # Adding a pseudocount of 1 to avoid issues with zeros

    # Only analyze common sample names between the two matrices
    clr_taxa = clr_taxa[, common_samples, drop = FALSE] # Subset clr_taxa to only include common samples
    clr_pathway = clr_pathway[, common_samples, drop = FALSE] # Subset clr_pathway to only include common samples

    common_samples2 = intersect(colnames(clr_pathway), colnames(clr_taxa)) # This wll print the common samples between the two matrices
    unique_samples2 = setdiff(colnames(clr_pathway), colnames(clr_taxa)) # This will print the samples that are unique to clr_pathway, in theory only week 7 should be different

    if (length(common_samples2) < 3) {
    stop("Not enough common samples (", length(common_samples2),
            ") between taxa and pathway CLR matrices for correlation. Exiting.")
    }
    else{
        print(paste0("Number of common samples between taxa and pathway CLR matrices for correlation: ", length(common_samples2)))
    }

    # Sanity checks for same sample order
    print("Looking at 20% random samples to check if the order of the samples in both matrices is the same")
    print("---- TAXA SAMPLE NAMES ----")
    random_samples = sample(1:ncol(clr_taxa), (length(1:ncol(clr_taxa)) / 5)) # Taking 20% random samples
    a = print(colnames(clr_taxa[, random_samples])) # Takes 20% random samples in the middle of the matrix
    print("---------------------------")
    print("---- PATHWAY SAMPLE NAMES ----")
    b = print(colnames(clr_pathway[, random_samples])) # Takes 20% random samples in the middle of the matrix
    print("---------------------------")

    if (all(a == b)) {
        print("The order of the samples in both matrices is the same, proceeding with correlation calculation: ")
    } else {
        print(".......WARNING.......")
        print("THE ORDER OF THE SAMPLES IN THE TWO MATRICES IS DIFFERENT, PLEASE CHECK.")
        print("EXITING FUNCTION")
        return(NULL)
    }

    # Here, I decide which self-correlation to do based on user input, taxa or pathways
    if(which_correlation == "taxa"){
        correlation_test = matrix(data = NA, nrow = nrow(clr_taxa), ncol = nrow(clr_taxa),
        dimnames = list(rownames(clr_taxa), rownames(clr_taxa))) # Setting dimension size and names to rows and columns
        p_values_mat = correlation_test
        for(i in seq_len(nrow(clr_taxa))){
            for(j in seq_len(nrow(clr_taxa))){
                test = cor.test(clr_taxa[i, ], clr_taxa[j, ], method = "spearman")
                correlation_test[i, j] = test$estimate
                p_values_mat[i, j]     = test$p.value
            }
        }
    }
    else{
        correlation_test = matrix(data = NA, nrow = nrow(clr_pathway), ncol = nrow(clr_pathway),
        dimnames = list(rownames(clr_pathway), rownames(clr_pathway))) # Setting dimension size and names to rows and columns
        p_values_mat = correlation_test
        for(i in seq_len(nrow(clr_pathway))){
            for(j in seq_len(nrow(clr_pathway))){
                test = cor.test(clr_pathway[i, ], clr_pathway[j, ], method = "spearman")
                correlation_test[i, j] = test$estimate
                p_values_mat[i, j]     = test$p.value
            }
        }
    }

    # P-value adjustment and removing NAs from rows
    pval_vector = p.adjust(as.vector(p_values_mat), method = "BH") # Convert matrix to vector and then do BH adjustment
    pval_matrix = matrix(pval_vector,
        nrow = nrow(p_values_mat), ncol = ncol(p_values_mat), # Convert back to matrix
        dimnames = list(rownames(p_values_mat), colnames(p_values_mat))
    )
    clean_correlation_test = correlation_test[rowSums(is.na(correlation_test)) != ncol(correlation_test), ] # Remove rows with all NAs for taxa
    clean_correlation_test = clean_correlation_test[, colSums(is.na(clean_correlation_test)) != nrow(clean_correlation_test)] # Remove columns with all NAs as well.

    pval_matrix_clean = pval_matrix[rownames(clean_correlation_test), colnames(clean_correlation_test)] # Subset pval_matrix to match clean_correlation_test

    # Subsetting to only include enriched pathways by Maaslin2
    significant_matrix = clean_correlation_test
    for (i in 1:nrow(significant_matrix)) {
        for (j in 1:ncol(significant_matrix)) {
            if (is.na(significant_matrix[i, j]) || abs(significant_matrix[i, j]) < rho_minimum || pval_matrix_clean[i, j] >= alpha_sign) {
                significant_matrix[i, j] = NA
            }
        }
    }

    sign_ass = sum(as.numeric(abs(clean_correlation_test) > rho_minimum & pval_matrix_clean < alpha_sign)) # This will give you the number of significant correlations with rho > 0.8 and adjusted p-value < 0.001
    range_ass = range(clean_correlation_test, na.rm = TRUE) # This will give you the range of correlation coefficients in your matrix
    
    print(paste0("Correlation test: ", dim(correlation_test)))
    print(paste0("Clean correlation test: ", dim(clean_correlation_test)))
    print(paste0("P-value matrix clean: ", dim(pval_matrix_clean)))
    print(paste0("Significant associations: ", sign_ass))
    print(paste0("Range of correlation coefficients: ", range_ass[1], " to ", range_ass[2]))

    # Printing analyses and files
    sentence4 = paste0("To make sure the order of the samples is the same, here are 6 random samples from each matrix:\n")
    sentence5 = paste0(a, " and ", b, "...............MATCH!\n")
    sentence6 = paste0(
        "A total of ", (nrow(correlation_test) - nrow(clean_correlation_test)), " taxa were removed  out of ", (nrow(correlation_test)), " due to having all NA Spearman correlation values, those are: ",
        paste(rownames(correlation_test)[which(rowSums(is.na(correlation_test)) == ncol(correlation_test))], collapse = ", "), "\n"
    )
    sentence7 = paste0(
        "A total of ", (ncol(correlation_test) - ncol(clean_correlation_test)), " pathways were removed  out of ", (ncol(correlation_test)), " due to having all NA Spearman correlation values, those are: ",
        paste(colnames(correlation_test)[which(colSums(is.na(correlation_test)) == nrow(correlation_test))], collapse = ", "), "\n"
    )
    sentence8 = paste0("There is a total of ", sign_ass, "/", length(clean_correlation_test), " significant associations with rho > ", rho_minimum, " and adjusted p-value < ", alpha_sign, "\n")
    sentence9 = paste0("The range of correlation coefficients is: ", range_ass[1], " to ", range_ass[2], "\n")
    writeLines(c(sentence4, sentence5, sentence6, sentence7, sentence8, sentence9), paste0(output_path, "/summary_analysis.txt"))
    write.table(clean_correlation_test, file = paste0(output_path, "/", filtering_value, "_", which_correlation, "_spearman_clean_correlation_test.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This saves a file where the columns/rows with NAs are removed from the Spearman correlation in taxa and pathways
    write.table(pval_matrix_clean, file = paste0(output_path, "/spearman_p_val_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This saves a file with the p-values instead of the spearman correlation rho after removing NAs from taxa and pathways
    write.table(significant_matrix, file = paste0(output_path, "/spearman_significant_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This prints the spearman values of those that have a greater 'n' rho, and lower 'n' q-value
    beepr::beep(1) # Beep when done
    
    return(clean_correlation_test)
}

normalizing_and_filtering = function(pathway_file, taxa_file, maaslin_filt_file, metadata_file, output_path, output_file_name_net, network_categories_output, 
filtering_category, filtering_value, alpha_sign, rho_minimum, ref_value) {
    # Open and read files 
    maaslin_filt = read.table(maaslin_filt_file, header = TRUE, sep = "\t")
    metadata_df = read.csv(metadata_file, header = TRUE, sep = ",", row.names = 1)
    pathway_mat = as.matrix(read.table(pathway_file, header = TRUE, sep = "\t", row.names = 1))
    taxa_mat = as.matrix(read.csv(taxa_file, header = TRUE, sep = ",", row.names = 1))
    # Renaming columns
    colnames(pathway_mat) = gsub("^X", "", colnames(pathway_mat)) # Removing any leading X from column names if present
    colnames(pathway_mat) = gsub("\\.", "-", colnames(pathway_mat)) # Replacing dots with dashes to match metacyc names
    rownames(taxa_mat) = gsub(" ", "_", rownames(taxa_mat)) # Replacing spaces with underscores in taxa names
    pathway_mat = t(pathway_mat)

    filtered_IDs = rownames(metadata_df[metadata_df[[filtering_category]] == filtering_value, , drop = FALSE])
    pathway_mat  = pathway_mat[, colnames(pathway_mat) %in% filtered_IDs, drop = FALSE]
    taxa_mat = taxa_mat[,colnames(taxa_mat) %in% filtered_IDs, drop = FALSE]

    pathway_mat = pathway_mat[rowSums(pathway_mat) > 0, ] # Remove pathways with all zeros
    taxa_mat = taxa_mat[rowSums(taxa_mat) > 0, ] # Remove taxa with all zeros

    common_samples = intersect(colnames(pathway_mat), colnames(taxa_mat)) # This wll print the common samples between the two matrices
    unique_samples = setdiff(colnames(pathway_mat), colnames(taxa_mat)) # This will print the samples that are unique to clr_pathway, in theory only week 7 should be different
    
    sample_sanity_check = sample_size_revision(common_samples, taxa_mat, pathway_mat)
    if(sample_sanity_check == TRUE){
        print("Proceeding with CLR normalization and correlation analysis...")
    }

    # CLR normalization
    clr_pathway = t(as.matrix(clr(t(pathway_mat) + 1))) # Adding a pseudocount of 1 to avoid issues with zeros
    clr_taxa = t(as.matrix(clr(t(taxa_mat) + 1))) # Adding a pseudocount of 1 to avoid issues with zeros

    # Only analyze common sample names between the two matrices
    clr_taxa = clr_taxa[, common_samples, drop = FALSE] # Subset clr_taxa to only include common samples
    clr_pathway = clr_pathway[, common_samples, drop = FALSE] # Subset clr_pathway to only include common samples

    common_samples2 = intersect(colnames(clr_pathway), colnames(clr_taxa)) # This wll print the common samples between the two matrices
    unique_samples2 = setdiff(colnames(clr_pathway), colnames(clr_taxa)) # This will print the samples that are unique to clr_pathway, in theory only week 7 should be different
    
    if (length(common_samples2) < 3) {
    stop("Not enough common samples (", length(common_samples2),
            ") between taxa and pathway CLR matrices for correlation. Exiting.")
    }
    else{
        print(paste0("Number of common samples between taxa and pathway CLR matrices for correlation: ", length(common_samples2)))
    }


    #clr_pathway = t(clr_pathway) # Transpose back to have pathways as rows and samples as columns
    #clr_taxa = t(clr_taxa) # Transpose back to have taxa as rows and samples as columns

    # Sanity checks for same sample order
    print("Looking at 20% random samples to check if the order of the samples in both matrices is the same")
    print("---- TAXA SAMPLE NAMES ----")
    random_samples = sample(1:ncol(clr_taxa), (length(1:ncol(clr_taxa)) / 5)) # Taking 20% random samples
    a = print(colnames(clr_taxa[, random_samples])) # Takes 20% random samples in the middle of the matrix
    print("---------------------------")
    print("---- PATHWAY SAMPLE NAMES ----")
    b = print(colnames(clr_pathway[, random_samples])) # Takes 20% random samples in the middle of the matrix
    print("---------------------------")

    if (all(a == b)) {
        print("The order of the samples in both matrices is the same, proceeding with correlation calculation: ")
    } else {
        print(".......WARNING.......")
        print("THE ORDER OF THE SAMPLES IN THE TWO MATRICES IS DIFFERENT, PLEASE CHECK.")
        print("EXITING FUNCTION")
        return(NULL)
    }

    correlation_test = matrix(              # Making a matrix with the same number of rows and columns for taxa and pathways. Each row is taxa, each column is a pathway
        data = NA, nrow = nrow(clr_taxa), ncol = nrow(clr_pathway),
        dimnames = list(rownames(clr_taxa), rownames(clr_pathway))
    )
    p_values_mat = correlation_test # Create an empty matrix for p-values with the same dimensions as correlation_test
    
    #Returns a vector of corrected names that are significant pathways under Maaslin2's eyes (q-val 0.001 and coef > 1 or < -1)
    corrected_pathways_names = correcting_names(maaslin_pathways = maaslin_filt, significant = "yes", filtering_category = filtering_category, filtering_value = filtering_value, ref_value = ref_value) #Modifying maaslin2 names for the right name and filtering by p-val and coefficient value Log2FC
    if(length(corrected_pathways_names) == 0){
        print("No significant pathways found by Maaslin2 under the given criteria, exiting function.")
        return(NULL)
    }
    print(colnames(correlation_test))
    print(corrected_pathways_names)
    # This part ensures that all the pathways found significant by Maaslin2 are present in the pathway matrix and that the naming is properly assigned in both matrices
    pathway_sanity_check = corrected_pathways_names %in% colnames(correlation_test)

    if(any(!pathway_sanity_check)){
        missing_pathways = corrected_pathways_names[!pathway_sanity_check] # Getting the names of the missing pathways 
        print(paste0("The following pathways were found significant by Maaslin2 but are missing in the pathway matrix: ", paste(missing_pathways, collapse = ", ")))
        stop("PATHWAY NAMING MISMATCH - EXITING FUNCTION")
    } 
    
    print(pathway_sanity_check)
    sentence1 = paste0("Originally, there were a total of ", nrow(pathway_mat), " pathways in the merged metagenome pathway prediction file from PICRUSt2 across all data sets. \n")
    sentence2 = paste0("After filtering by significant pathways found in the ", basename(taxa_file), " community, there are a total of: ", length(corrected_pathways_names), " pathways. \n")
    sentence3 = paste0("The number of pathways calculated by Maaslin2 in the ", maaslin_filt_file, " is ", nrow(maaslin_filt))
    print(sentence1)
    print(sentence2)
    print(sentence3)

    # Running the spearman correlation test and storing values in the empty matrices
    for (i in seq_len(nrow(clr_taxa))) {        # over taxa
        for (j in seq_len(nrow(clr_pathway))) {   # over pathways
        test = cor.test(clr_taxa[i, ], clr_pathway[j, ], method = "spearman")
        correlation_test[i, j] = test$estimate
        p_values_mat[i, j]     = test$p.value
        }
    }

    # P-value adjustment and removing NAs from rows
    pval_vector = p.adjust(as.vector(p_values_mat), method = "BH") # Convert matrix to vector and then do BH adjustment
    pval_matrix = matrix(pval_vector,
        nrow = nrow(p_values_mat), ncol = ncol(p_values_mat), # Convert back to matrix
        dimnames = list(rownames(p_values_mat), colnames(p_values_mat))
    )
    clean_correlation_test = correlation_test[rowSums(is.na(correlation_test)) != ncol(correlation_test), ] # Remove rows with all NAs for taxa
    clean_correlation_test = clean_correlation_test[, colSums(is.na(clean_correlation_test)) != nrow(clean_correlation_test)] # Remove columns with all NAs for pathways

    # Subsetting to only include enriched pathways by Maaslin2
    enriched_clean_correlation_test = clean_correlation_test[, colnames(clean_correlation_test) %in% corrected_pathways_names]
    pval_matrix_clean = pval_matrix[rownames(clean_correlation_test), colnames(clean_correlation_test)] # Subset pval_matrix to match clean_correlation_test

    print(paste0("There are ", paste(dim(clean_correlation_test), collapse = " and "), " taxa and pathways in the rho correlation matrix, respectively."))
    print(paste0("There are ", paste(dim(pval_matrix_clean), collapse = " and "), " taxa and pathways in the clean p-value matrix, respectively."))
    print(paste0("There are ", paste(dim(enriched_clean_correlation_test), collapse = " and "), " taxa and enriched pathways in the clean correlation test, respectively."))

    # Filtering by significance and rho value
    significant_matrix = clean_correlation_test
    for (i in 1:nrow(significant_matrix)) {
        for (j in 1:ncol(significant_matrix)) {
            if (is.na(significant_matrix[i, j]) || abs(significant_matrix[i, j]) < rho_minimum || pval_matrix_clean[i, j] >= alpha_sign) {
                significant_matrix[i, j] = NA
            }
        }
    }

    sign_ass = sum(as.numeric(abs(clean_correlation_test) > rho_minimum & pval_matrix_clean < alpha_sign)) # This will give you the number of significant correlations with rho > 0.8 and adjusted p-value < 0.001
    range_ass = range(clean_correlation_test, na.rm = TRUE) # This will give you the range of correlation coefficients in your matrix


    # Printing analyses and files
    sentence4 = paste0("To make sure the order of the samples is the same, here are 6 random samples from each matrix:\n")
    sentence5 = paste0(a, " and ", b, "...............MATCH!\n")
    sentence6 = paste0(
        "A total of ", (nrow(correlation_test) - nrow(clean_correlation_test)), " taxa were removed  out of ", (nrow(correlation_test)), " due to having all NA Spearman correlation values, those are: ",
        paste(rownames(correlation_test)[which(rowSums(is.na(correlation_test)) == ncol(correlation_test))], collapse = ", "), "\n"
    )
    sentence7 = paste0(
        "A total of ", (ncol(correlation_test) - ncol(clean_correlation_test)), " pathways were removed  out of ", (ncol(correlation_test)), " due to having all NA Spearman correlation values, those are: ",
        paste(colnames(correlation_test)[which(colSums(is.na(correlation_test)) == nrow(correlation_test))], collapse = ", "), "\n"
    )
    sentence8 = paste0("There is a total of ", sign_ass, "/", length(clean_correlation_test), " significant associations with rho > ", rho_minimum, " and adjusted p-value < ", alpha_sign, "\n")
    sentence9 = paste0("The range of correlation coefficients is: ", range_ass[1], " to ", range_ass[2], "\n")
    writeLines(c(sentence1, sentence2, sentence3, sentence4, sentence5, sentence6, sentence7, sentence8, sentence9), paste0(output_file_name_net, "_summary_analysis.txt"))

    write.table(clean_correlation_test, file = paste0(output_file_name_net,"_spearman_clean_correlation_test.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This saves a file where the columns/rows with NAs are removed from the Spearman correlation in taxa and pathways
    write.table(pval_matrix_clean, file = paste0(output_file_name_net, "_spearman_p_val_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This saves a file with the p-values instead of the spearman correlation rho after removing NAs from taxa and pathways
    write.table(significant_matrix, file = paste0(output_file_name_net, "_spearman_significant_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This prints the spearman values of those that have a greater 'n' rho, and lower 'n' q-value
    write.table(enriched_clean_correlation_test, file = paste0(output_file_name_net, "_spearman_maaslin2_enriched_clean_correlation_test.tsv"), sep = "\t", quote = FALSE, col.names = NA) # This one writes a file with the spearman correlation only found in the pathways enriched in maaslin2 when comparing it against another group.
    
    writeLines(colnames(enriched_clean_correlation_test), paste0(network_categories_output,".txt")) # This writes a file with the names of the pathways that were found significant by Maaslin2 and used for the network construction
    
    beepr::beep(1) # Beep when done
    return(enriched_clean_correlation_test)
}

make_encoded_matrix = function(per_seq_contrib, asv_to_taxa_df, taxa_list, pathway_list, min_copy = 1e-9){
    # Read the collapsed contribution file (from mapping_taxa_to_pathway function)
    df_seq_contribution = read.table(per_seq_contrib, header = TRUE, sep = "\t")
    
    # Filter to only include taxa and pathways from the correlation analysis
    df_seq_contribution = df_seq_contribution[df_seq_contribution$genus_and_species %in% taxa_list & df_seq_contribution$function. %in% pathway_list, ] 
    
    # Determine if each species encodes the pathway (genome_function_count > min_copy)
    df_seq_contribution$encodes = df_seq_contribution$genome_function_count > min_copy 
    
    # Get unique pairs of pathways and taxa where encoding occurs
    pairs_correlated = unique(df_seq_contribution[df_seq_contribution$encodes, c("function.", "genus_and_species")]) 

    # Create empty matrix
    encoded_matrix = matrix(FALSE, nrow = length(taxa_list), 
                                ncol = length(pathway_list),
                                dimnames = list(taxa_list, pathway_list))
    
    # Fill matrix with TRUE where species encodes pathway
    if (nrow(pairs_correlated) > 0){
        encoded_matrix[cbind(match(pairs_correlated$genus_and_species, taxa_list),
                            match(pairs_correlated$function., pathway_list))] <- TRUE
    } 
    View(encoded_matrix)
    return(encoded_matrix)
}

heatmap_plotting = function(input_data_name, matrix_to_plot, per_sequence_contribution_file){
    if(length(matrix_to_plot) == 0){
        print("No data available for heatmap plotting, exiting function.")
        return(NULL)
    }
    matrix_to_plot <- matrix_to_plot[!grepl("_NA", rownames(matrix_to_plot)), , drop = FALSE] # Remove rows with '_NA' in their names
    matrix_to_plot <- matrix_to_plot[, !grepl("_NA", colnames(matrix_to_plot)), drop = FALSE] # Remove columns with '_NA' in their names
    markdown = "X"
    encoded_mat = make_encoded_matrix(per_seq_contrib = per_sequence_contribution_file,
        asv_to_taxa_df = NULL,
        taxa_list = rownames(matrix_to_plot),
        pathway_list = colnames(matrix_to_plot),
        min_copy = 1e-9
    )

    png(paste0(input_data_name, ".png"), width = 2800, height = 2400, res = 150)
    my_palette = colorRampPalette(c("purple", "blue", "white", "orange", "darkred"))(n = 30)
    
    # Create cellnote matrix with black circle where encoding is TRUE
    cellnote_matrix = ifelse(encoded_mat, "●", "")

    heatmap1 = heatmap.2(matrix_to_plot,
        col = my_palette,
        trace = "none",
        density.info = "none",
        dendrogram = "both",
        cellnote = cellnote_matrix,  # Add circles for encoded pathways
        notecex = 2.5,  # Size of circles
        notecol = "black",  # Color of circles
        key.title = "Spearman\nCorrelation",
        key.xlab = "Coefficient (ρ)",
        key.par = list(cex = 1.1), # Legend text size
        keysize = 0.8, # Legend size
        margins = c(40, 20), # Margins
        cexRow = 2.1, # Row text
        cexCol = 1.3, # Column text
        colsep = 0,
        rowsep = 0:nrow(matrix_to_plot),
        sepcolor = "black",
        lmat = rbind(c(0, 3, 4), c(2, 1, 0)),
        sepwidth = c(0.04, 0.04),
        lhei = c(0.6, 2.5), lwid = c(0.6, 4, 1.5), lwd = 2
    )
    taxa_ordered = rownames(matrix_to_plot)[heatmap1$rowInd] # Print row names in the order they appear in the heatmap
    pathways_ordered = colnames(matrix_to_plot)[heatmap1$colInd] # Print column names in the order they appear in the heatmap
    encoded_ordered = encoded_mat[taxa_ordered, pathways_ordered, drop = FALSE] # Reorder the encoded matrix to match heatmap order
    
    not_encoded = which(encoded_ordered == FALSE, arr.ind = TRUE)
    if(nrow(not_encoded) > 0){
        x = not_encoded[, "col"]
        y = not_encoded[, "row"]
        text(x,y, labels = markdown, cex = 0.9, font = 2, col = "black")
    }

    write.table(colnames(matrix_to_plot), file = paste0(input_data_name,".tsv"), sep = "\t", quote = FALSE, col.names = FALSE)

    matrix_col_names_sorted = colnames(matrix_to_plot)[heatmap1$colInd] # Print column names in the order they appear in the heatmap
    write.table(matrix_col_names_sorted, file = paste0(input_data_name,"_sorted.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
    
    dev.off()
    beep(1)
    return(NULL)
}

retrieving_correlated_features = function(){
}

mapping_taxa_to_pathway = function(picrust2_file_strat, asv_to_taxa, consortia){
    df_picrust2 = read.table(picrust2_file_strat, header = TRUE, sep = "\t")
    df_asv_to_taxa = read.table(asv_to_taxa, header = TRUE, sep = "\t", row.names = 1)
    df_asv_to_taxa$ASV_ID = rownames(df_asv_to_taxa) # Moving rownames to a column
    df_asv_to_taxa = df_asv_to_taxa[, c("ASV_ID", setdiff(colnames(df_asv_to_taxa), "ASV_ID"))] # Putting the ASV_ID column first
    df_picrust2 = df_picrust2 %>% rename("ASV_ID" = "sequence")
    df_merged = full_join(df_picrust2, df_asv_to_taxa, by = "ASV_ID")
    df_merged$genus_and_species = paste0(df_merged$genus_final, "_", df_merged$species_final)
    df_merged = na.omit(df_merged)
    
    # Collapse ASVs from the same species contributing to the same pathway by summing their contributions
    df_collapsed = df_merged %>%
        group_by(sample, function., genus_and_species) %>%
        summarise(
            taxon_abun = sum(taxon_abun, na.rm = TRUE),
            taxon_rel_abun = sum(taxon_rel_abun, na.rm = TRUE),
            genome_function_count = sum(genome_function_count, na.rm = TRUE),
            taxon_function_abun = sum(taxon_function_abun, na.rm = TRUE),
            taxon_rel_function_abun = sum(taxon_rel_function_abun, na.rm = TRUE),
            norm_taxon_function_contrib = sum(norm_taxon_function_contrib, na.rm = TRUE),
            n_ASVs = n(),
            .groups = 'drop'
        )
    
    write.table(df_merged, file = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/",consortia,"_inocula_output/",consortia,"_pathway_out_contrib/",consortia,"_picrust2_pathway_with_taxa_names_uncollapsed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(df_collapsed, file = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/",consortia,"_inocula_output/",consortia,"_pathway_out_contrib/",consortia,"_picrust2_pathway_with_taxa_names_collapsed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    return(df_collapsed)
}

self_network = function(pathway_file_name, taxa_file_name, metadata_file_name, output_path, rho_minimum, filtering_category, filtering_value,
                        community, timepoint, clean_correlation_matrix_file, output_file){
# Connecting to Cytoscape
print(cytoscapePing()) # Checking connection to Cytoscape

# Loading files
    taxa_file = read.csv(taxa_file_name, header = TRUE, sep = ",", row.names = 1)
    path_file = as.matrix(read.table(pathway_file_name, header = TRUE, sep = "\t", row.names = 1))
    colnames(path_file) = gsub("^X", "", colnames(path_file)) # Removing any leading X from column names if present
    colnames(path_file) = gsub("\\.", "-", colnames(path_file)) # Replacing dots with dashes to match metacyc names
    rownames(taxa_file) = gsub(" ", "_", rownames(taxa_file)) # Replacing spaces with underscores in taxa names
    path_file = t(path_file)
    metadata_df = read.csv(metadata_file_name, header = TRUE, sep = ",", row.names = 1)
    taxa_file = t(taxa_file)

    # Skipping over the Maaslin2 analysis, as this is a self-correlation network

    metadata_filtered = metadata_df[(metadata_df$Merged_weeks) == timepoint, ]
    metadata_filtered = metadata_filtered[metadata_filtered$Subcommunity == community, ]
    
    path_matrix_clean = path_file[,colnames(path_file) %in% metadata_filtered$ID, drop = FALSE] #Filtering pathway matrix by IDs matching the community and timepoint of the filtered metadata
    taxa_matrix_clean = as.matrix(taxa_file[,colnames(taxa_file) %in% metadata_filtered$ID, drop = FALSE]) #Filtering by sample IDs matching the community and timepoint
    total_mice = ncol(path_matrix_clean) #Getting total number of mice after filtering by metadata filtered variable, 52


}

make_network = function(pathway_file_name, taxa_file_name, maaslin_file_name, metadata_file_name, output_path, rho_minimum, filtering_category, filtering_value,
                        community, timepoint, clean_correlation_matrix_file, ref_value, output_file) {
    #Connecting to Cytoscape
    print(cytoscapePing()) #Checking connection to Cytoscape

    #Loading files
    taxa_file = read.csv(taxa_file_name, header = TRUE, sep = ",", row.names = 1)
    path_file = as.matrix(read.table(pathway_file_name, header = TRUE, sep = "\t", row.names = 1))
    colnames(path_file) = gsub("^X", "", colnames(path_file)) # Removing any leading X from column names if present
    colnames(path_file) = gsub("\\.", "-", colnames(path_file)) # Replacing dots with dashes to match metacyc names
    rownames(taxa_file) = gsub(" ", "_", rownames(taxa_file)) # Replacing spaces with underscores in taxa names
    path_file = t(path_file)
    maaslin_filt = read.table(maaslin_file_name, header = TRUE, sep = "\t")
    metadata_df = read.csv(metadata_file_name, header = TRUE, sep = ",", row.names = 1)
    taxa_file = t(taxa_file)

    # Analyzing and retrieving significant pathways from Maaslin2 and names corrected for the pathway matrix
    corrected_names = correcting_names(maaslin_pathways = maaslin_filt, significant = "yes", filtering_category = filtering_category, filtering_value = filtering_value, ref_value = ref_value) #This fixes the naming of pathways in Maaslin2 from the file in question (NS1 vs S2 w9w10), and takes all the pathways (significant or not)
    if(length(corrected_names) == 0){
        print("No significant pathways found by Maaslin2 under the given criteria, exiting function.")
        return(NULL)
    }
    path_matrix_significant = path_file[rownames(path_file) %in% corrected_names, ] #This filters the WHOLE pathway matrix of all samples and filters by the ones found in the Maaslin2 comparison

    #Filtering by consortia and by timepoint:
    metadata_filtered = metadata_df[metadata_df[[filtering_category]] == filtering_value, , drop = FALSE]

    path_matrix_clean = path_matrix_significant[,colnames(path_matrix_significant) %in% metadata_filtered$ID, drop = FALSE] #Filtering pathway matrix by IDs matching the community and timepoint of the filtered metadata
    taxa_matrix_clean = as.matrix(taxa_file[,colnames(taxa_file) %in% metadata_filtered$ID, drop = FALSE]) #Filtering by sample IDs matching the community and timepoint
    total_mice = ncol(path_matrix_clean) #Getting total number of mice after filtering by metadata filtered variable, 52

    pathway_without_mouse = nrow(path_matrix_clean[rowSums(path_matrix_clean) == 0, ]) #Identifies number of pathways that are zero in all mice, 4
    taxa_without_mouse = nrow(taxa_matrix_clean[rowSums(taxa_matrix_clean) == 0, ]) #Identifies number of taxa that are zero in all mice, 6

    print(paste0("Pathways removed due to zero abundance: ", paste(rownames(path_matrix_clean)[rowSums(path_matrix_clean) == 0], collapse = ", ")))

    #These two will be used for the network and pivotted tables
    path_matrix_final = path_matrix_clean[rowSums(path_matrix_clean) > 0, , drop = FALSE] # This command filters the pathway matrix to only include those pathways that have at least one mouse present
    taxa_matrix_final = taxa_matrix_clean[rowSums(taxa_matrix_clean) > 0, , drop = FALSE] #31 rows are taxa, 52 columns are mice

    num_of_mice_in_taxa = ((rowSums(taxa_matrix_final > 0))/(ncol(taxa_matrix_final)))*100 #Getting the proportional number of mice present in each taxa
    num_of_mice_in_pathway = ((rowSums(path_matrix_final > 0))/(ncol(path_matrix_final)))*100 #Getting the proportional number of mice present in each pathway

    log_file = file.path(paste0(output_file,"_log.txt"))
    con = file(log_file, open = "wt")

    sentence1 = "NETWORK ANALYSIS SUMMARY: "
    sentence2 = "--------------------------------------------------"
    sentence3 = paste0("The original pathway file with all samples and pathways has ", nrow(path_file), " pathways and ", ncol(path_file), " mice samples.")
    sentence4 = paste0("After filtering by significant pathways using the Maaslin2 output file ", maaslin_file_name, ", there are a total of ", nrow(path_matrix_significant), " pathways in ", ncol(path_matrix_significant), " mice samples.")
    sentence5 = paste0("Finally, after filtering the pathway file by only including those with at least one mouse present, and by metadata chosen categories, there are a total of ", nrow(path_matrix_final), " pathways in ", ncol(path_matrix_final), " mice samples.\n")
    sentence6 = paste0("The original taxa file for the analyzed community has ", nrow(taxa_file), " taxa and ", ncol(taxa_file), " mice samples in ", community, ".")
    sentence7 = paste0("Finally, after filtering the taxa file by only including those with at least one mouse present, and by metadata chosen categories, there are a total of ", ncol(taxa_matrix_final), " taxa in ", nrow(taxa_matrix_final), " mice samples.\n")
    sentence8 = paste0("TOTAL NUMBER OF FINAL TAXA TO ANALYZE: ", nrow(taxa_matrix_final), "\n TOTAL NUMBER OF FINAL PATHWAYS TO ANALYZE: ", nrow(path_matrix_final), "\n TOTAL NUMBER OF MICE SAMPLES TO ANALYZE: ", total_mice, "\n")
    sentence9 = ("--------------------------------------------------")

    sentences = c(sentence1, sentence2, sentence3, sentence4, sentence5, sentence6, sentence7, sentence8, sentence9)
    writeLines(sentences, con = con)

    df_rho_estimate = as.data.frame(clean_correlation_matrix_file) #Using the cleaned correlation matrix from the previous function
    df_rho_estimate = rownames_to_column(df_rho_estimate, var = "taxon") #This converts the rownames (taxa) into a column called 'taxon' in the dataframe df_rho_estimate

    pivoted_rho_table = pivot_longer(data = df_rho_estimate, cols = -taxon,
                                names_to = "pathway", values_to = "rho") #This 'pivots' (converts from wide to long format) the correlation matrix by taking all columns and converting them into two columns: 'pathway' and 'rho' (correlation value)

    pivoted_edge_table = pivoted_rho_table %>%
        mutate(filtered_rho = ifelse((abs(rho)) < rho_minimum, 0, (rho))) %>% #Filtering by rho value. If its absolute value is less than the minimum (my case, 0.70), set it to 0
        mutate(sign_of_rho = ifelse(rho > 0, "positive", "negative"), corr_value = abs(rho)) %>%
        mutate(source = gsub("_", " ", taxon), target = gsub("_", " ", pathway), interaction = sign_of_rho, weight = rho) %>% #Adding sign of rho and converting to absolute value for node size later
        mutate(taxon_and_pathway = paste0(gsub("_", " ", taxon), "&", gsub("_", " ", pathway))) #Creating a new column with the taxon and pathway names combined for easier filtering later

    pivoted_edge_table = pivoted_edge_table %>%
        mutate(source = gsub("_", " ", source)) %>%
        mutate(taxon_and_pathway = gsub("_", " ", taxon_and_pathway))
    
    #Sanity check(looking at a random rows), choosing 10 random numbers from a uniform distribution to test if the values match in the pivoted table and the original matrices.
    for(i in round(runif(10, min = 1, max = nrow(pivoted_edge_table)))) {
        print(paste0("Row ", i, ": Taxon = ", pivoted_edge_table$taxon[i], ", Pathway = ", pivoted_edge_table$pathway[i], ", Rho = ", pivoted_edge_table$rho[i]))
        rho_original_value = df_rho_estimate[df_rho_estimate$taxon == pivoted_edge_table$taxon[i], pivoted_edge_table$pathway[i]]
        if(rho_original_value == pivoted_edge_table$weight[i]) {
            writeLines(paste0("Rho matches original matrix: ", rho_original_value, " == ", pivoted_edge_table$rho[i]), con = con)
        } else {
            writeLines("RHO VALUE DOES NOT MATCH ORIGINAL MATRIX!", con = con)
            return(NULL)
        }
        writeLines("--------------------------------------------------", con = con)
    }

    degree_taxa = pivoted_edge_table %>% dplyr::count(taxon, wt = (filtered_rho != 0), name = "degree") #Centrality degree calculation in taxa
    degree_pathway = pivoted_edge_table %>% dplyr::count(pathway, wt = (filtered_rho != 0), name = "degree") #Centrality degree calculation in pathways

    degree_taxa = degree_taxa %>% mutate(taxon = gsub("_", " ", taxon)) #Replacing underscores with spaces for better visualization in Cytoscape
    degree_pathway = degree_pathway %>% mutate(pathway = gsub("_", " ", pathway)) #Replacing underscores with spaces for better visualization in Cytoscape
    
    degree_node = bind_rows(degree_taxa, degree_pathway) #Merging both into one

    pivoted_edge_table = pivoted_edge_table %>% #Filters the edge table to only include those with filtered rho not equal to 0, meaning those that passed the rho minimum threshold
        filter(filtered_rho != 0)

    # Identify pathways excluded due to weak correlations
    all_maaslin_pathways = gsub("_", " ", colnames(clean_correlation_matrix_file))
    pathways_in_network = unique(pivoted_edge_table$target)
    excluded_pathways = setdiff(all_maaslin_pathways, pathways_in_network)
    
    if(length(excluded_pathways) > 0) {
        writeLines(paste0("\nPathways excluded from network due to no correlations meeting rho >= ", rho_minimum, " threshold:"), con = con)
        writeLines(paste0("Total excluded: ", length(excluded_pathways), " out of ", length(all_maaslin_pathways), " Maaslin2-significant pathways"), con = con)
        writeLines(paste0("Excluded pathways: ", paste(excluded_pathways, collapse = ", ")), con = con)
        writeLines("--------------------------------------------------", con = con)
    }

    nodes_ids_in_both = sort(unique(c(pivoted_edge_table$source, pivoted_edge_table$target))) #Getting unique node IDs present in the edge table after filtering by rho minimum
    close(con)

    taxa_ids     = unique(pivoted_edge_table$source)   
    pathway_ids  = unique(pivoted_edge_table$target) 

    nodes_ids_in_both = sort(unique(c(pivoted_edge_table$source, pivoted_edge_table$target)))

    nodes_table = tibble(id = nodes_ids_in_both) %>%
                                mutate(
                                    group = case_when(
                                        id %in% taxa_ids ~ "taxa",
                                        id %in% pathway_ids ~ "pathway",
                                        TRUE ~ "unknown"
                                    ),
                                    score = case_when(
                                        group == "taxa" ~ {
                                            id_with_underscore = gsub(" ", "_", id)
                                            as.numeric(num_of_mice_in_taxa[id_with_underscore])
                                        },
                                        group == "pathway" ~ {
                                            id_with_underscore = gsub(" ", "_", id)
                                            as.numeric(num_of_mice_in_pathway[id_with_underscore])
                                        },
                                        TRUE ~ NA_real_
                                    )
                                ) %>%
                                left_join(degree_taxa, by = c("id" = "taxon")) %>%
                                left_join(degree_pathway, by = c("id" = "pathway")) %>%
                                mutate(degree = coalesce(degree.x, degree.y)) %>%
                                select(-degree.x, -degree.y) %>%
                                mutate(id = gsub("_", " ", id))
    
    nodes_table = nodes_table %>% mutate(id = gsub("_", " ", id)) #Replacing underscores with spaces for better visualization in Cytoscape

    nodes_table = nodes_table %>%
    mutate(opacity_category = case_when(
        score >= 0 & score < 20 ~ "1-20%",
        score >= 20 & score < 40 ~ "20-40%", 
        score >= 40 & score < 60 ~ "40-60%",
        score >= 60 & score < 80 ~ "60-80%",
        score >= 80 & score <= 100 ~ "80-100%",
        TRUE ~ "Unknown"
    ),
    size_category = case_when(
        degree >= 1 & degree < 3 ~ "1-2",
        degree >= 3 & degree < 5 ~ "3-4", 
        degree >= 5 & degree < 7 ~ "5-6",
        degree >= 7 & degree < 8 ~ "7-8",
        degree >= 8 & degree < 10000 ~ "8+",
        TRUE ~ "Unknown"
    ))

    #Generating the network in Cytoscape
    print("---------------------------------------------------")
    print(paste0("Creating network ", output_file, " in Cytoscape. Please be patient."))

    name_for_network = substr(basename(output_file), 1, nchar(basename(output_file))) #Getting just the last part of the path and removing the _log.txt extension
    createNetworkFromDataFrames(nodes_table, pivoted_edge_table,
                                title = name_for_network)

    #Setting default visual styles for the network
    #Nodes first
    setNodeShapeDefault("ellipse", style.name = "default")             
    setNodeFontSizeDefault(new.size = 18, style.name = "default")       #Setting default font size to 18
    setNodeLabelPositionDefault(new.nodeAnchor = "S",                   #Setting position of the node and text label
                new.graphicAnchor = "S", new.justification = "c", 
                new.xOffset = 0, new.yOffset = 20, style.name = "default")
    setNodeBorderWidthDefault(new.width = 2, style.name = "default")    # Setting border width to 2 for all nodes
    setNodeBorderColorDefault(new.color = "black", style.name = "default")  #Setting border color to black for all nodes
    setNodeLabelColorDefault(new.color = "white", style.name = "default")                          #Setting default font color to white
    setNodeFillOpacityMapping(table.column = "opacity_category",
                        table.column.values = c("1-20%", "20-40%", "40-60%", "60-80%", "80-100%"),
                        opacities = c(25, 45, 85, 125, 225),
                        mapping.type = "discrete",  # Changed to discrete for categories
                        style.name = "default")
    setNodeColorMapping(table.column = "group",     #Taxa will have green nodes, pathways orange.
                    table.column.values = c("taxa","pathway"),
                    colors = c("darkblue", "#ca7717"),
                    mapping.type = "discrete",
                    style.name = "default")
    
    setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_BACKGROUND_SHAPE", 
                            value = "rectangle"), style.name = "default")   # Shape of the label_background
    setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_BACKGROUND_COLOR", 
                            value = "black"), style.name = "default")       # label background to black
    setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_BACKGROUND_OPACITY", 
                            value = 250), style.name = "default")           # label background opacity to none
    

    longnames = any(nchar(nodes_table$id) > 19) #Check if any node name is longer than 15 characters
    
    if(longnames == TRUE){
       setNodeLabelPositionDefault(new.nodeAnchor = "S",
                               new.graphicAnchor = "N", 
                               new.justification = "c",
                               new.xOffset = 0, 
                               new.yOffset = 10,  # Move further down
                               style.name = "default")
    } else {
    # For normal names, use standard position
    setNodeLabelPositionDefault(new.nodeAnchor = "S",
                               new.graphicAnchor = "N", 
                               new.justification = "c",
                               new.xOffset = 0, 
                               new.yOffset = 0,  # Standard position
                               style.name = "default")
    }
    setNodeSizeMapping(table.column = "size_category", #The degree of the node defines its size, in this case, the number of connections in a node
                    table.column.values = c("1-2", "3-4", "5-6", "7-8", "8+"),
                    sizes = c(20, 40, 60, 80, 100),
                    mapping.type = "discrete",
                    style.name = "default")
    setEdgeLineWidthMapping(table.column = "corr_value", #The correlation value is the absolute number from rho, defines the thickness of the edge by rho value
                        table.column.values = c((rho_minimum), max(pivoted_edge_table$corr_value)),
                        widths = c(0.5, 4),
                        mapping.type = "continuous",
                        style.name = "default")
    setEdgeColorMapping(table.column = "interaction", #If rho is positive, then the edge will be black. If negative, then red.
                        table.column.values= c("positive", "negative"),
                        colors = c("black", "red"),
                        mapping.type = "discrete",
                        style.name = "default")

    beep(8) #Beep when done

}

global = function() {
    loading_packages()

    decision = 2 #1 for single community analysis, 2 for pairwise community analysis, 3 for self-comparison
    inocula = TRUE

    consortia = c("NS1","S2")
    timepoints = c("w5", "w9w10")
    pathway_file = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/Pathway_merged_metagenome.tsv"
    metadata_file = "/Users/danielcm/Desktop/diammatics/T1D/metadata_ps_final.csv"
    filtering_category = "Week_and_consortia"
    alpha_sign = 0.05
    rho_minimum = 0.7
    heatmap_output_path = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Heatmaps/"
    network_output_path = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/"
    picrust2_path = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/"

    if(decision == 1){
        for(i in seq(consortia)){
            for(j in seq(timepoints)){
                taxa_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ps_",tolower(consortia[[i]]),"_final.csv")
                maaslin_filt_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/pairwise_comparisons/same_community/Pathway/maaslin2_Week_and_consortia_pathway_", consortia[[i]], "_", timepoints[[1]], "_vs_", 
                                        consortia[[i]], "_", timepoints[[2]], "_ref_Week_and_consortia,", consortia[[i]], "_", timepoints[[2]], "/all_results.tsv")
                filtering_value = paste0(consortia[[i]], "_", timepoints[[j]])
                output_path = paste0("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/", consortia[[i]])
                output_file_name_heat = paste0(heatmap_output_path, consortia[[i]],"_", timepoints[[j]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[1]], "_and_",consortia[[i]], "_", timepoints[[2]])
                output_file_name_net = paste0(output_path, "/", consortia[[i]], "_", timepoints[[j]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[1]], "_and_",consortia[[i]],"_", timepoints[[2]])
                output_file_ids_network = paste0(output_path, "/", consortia[[i]], "_", timepoints[[j]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[1]], "_and_",consortia[[i]],"_", timepoints[[2]], "_network_pathway_IDs")
                ref_value = sub(".*,(.*?)\\/.*", "\\1", maaslin_filt_file) # Extracting the reference value from the maaslin filtered file path)
                
                enriched_clean_correlation_test = normalizing_and_filtering(
                    pathway_file = pathway_file,
                    taxa_file = taxa_file,
                    maaslin_filt_file = maaslin_filt_file,
                    metadata_file = metadata_file,
                    output_path = output_path,
                    filtering_category = filtering_category,
                    filtering_value = filtering_value,
                    alpha_sign = as.numeric(alpha_sign),
                    rho_minimum = as.numeric(rho_minimum),
                    output_file_name_net = output_file_name_net,
                    network_categories_output = output_file_ids_network,
                    ref_value = ref_value)

                #heatmap_plotting(input_data_name = output_file_name_heat, matrix_to_plot = enriched_clean_correlation_test, inoc1, inoc2)

                make_network(pathway_file_name = pathway_file,
                    taxa_file_name = taxa_file,
                    maaslin_file_name = maaslin_filt_file,
                    metadata_file_name = metadata_file,
                    output_path = output_path,
                    filtering_category = filtering_category,
                    filtering_value = filtering_value,
                    rho_minimum = as.numeric(rho_minimum),
                    community = consortia[[i]],
                    timepoint = timepoints[[j]],
                    clean_correlation_matrix_file = enriched_clean_correlation_test,
                    ref_value = ref_value,
                    output_file = output_file_name_net)
            }
        }
    }

    else if(decision == 2){
        m_values = c(1, 2) #This variable will help alternate between the two consortia in the pairwise comparison as reference and filtering value
        for(i in 1:length(consortia)){
            for(j in 1:length(consortia)){
                if(i==j | j<i){
                    next
                }   
                for(k in 1:length(timepoints)){
                    for(m in seq(m_values)){
                        print("-----------------------------------")
                        print(paste0("PROCESSING CONSORTIA ", consortia[[i]], " AND ", consortia[[j]], " AT TIMEPOINT ", timepoints[[k]], " WITH REFERENCE AS ", ifelse(m == 1, consortia[[i]], consortia[[j]])))

                        maaslin_filt_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/pairwise_comparisons/",timepoints[[k]],"/Pathway/maaslin2_Week_and_consortia_pathway_", consortia[[i]], "_", timepoints[[k]],"_vs_", 
                                                consortia[[j]], "_", timepoints[[k]], "_ref_Week_and_consortia,", consortia[[j]], "_", timepoints[[k]], "/all_results.tsv")
                        if(m == 1){
                            filtering_value = paste0(consortia[[i]],"_", timepoints[[k]])
                            output_file_name_heat = paste0(heatmap_output_path,"/", consortia[[i]],"_", timepoints[[k]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[k]], "_and_",consortia[[j]], "_", timepoints[[k]])
                            output_path = paste0("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/", consortia[[i]])
                            output_file_name_net = paste0(output_path, "/", consortia[[i]], "_", timepoints[[k]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[k]], "_and_",consortia[[j]],"_", timepoints[[k]])
                            taxa_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ps_",tolower(consortia[[i]]),"_final.csv")
                            picrust2_file_strat = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/",consortia[[i]],"_inocula_output/",consortia[[i]],"_pathway_out_contrib/path_abun_contrib.tsv")
                            asv_to_taxa = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ASV_to_taxa_",consortia[[i]],"_inoc.tsv")
                            community = consortia[[i]]
                        }
                        else if(m == 2){
                            filtering_value = paste0(consortia[[j]],"_", timepoints[[k]])
                            output_path = paste0("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/", consortia[[j]])
                            output_file_name_heat = paste0(heatmap_output_path, consortia[[j]],"_", timepoints[[k]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[k]], "_and_",consortia[[j]], "_", timepoints[[k]])
                            output_file_name_net = paste0(output_path, "/", consortia[[j]], "_", timepoints[[k]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[k]], "_and_",consortia[[j]],"_", timepoints[[k]])
                            taxa_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ps_",tolower(consortia[[j]]),"_final.csv")
                            picrust2_file_strat = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/",consortia[[j]],"_inocula_output/",consortia[[j]],"_pathway_out_contrib/path_abun_contrib.tsv")
                            asv_to_taxa = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ASV_to_taxa_",consortia[[j]],"_inoc.tsv")
                            community = consortia[[j]]
                        }
                        ref_value = sub(".*,(.*?)\\/.*", "\\1", maaslin_filt_file) # Extracting the reference value from the maaslin filtered file path)
                        
                        enriched_clean_correlation_test = normalizing_and_filtering(
                            pathway_file = pathway_file,
                            taxa_file = taxa_file,
                            maaslin_filt_file = maaslin_filt_file,
                            metadata_file = metadata_file,
                            output_path = output_path,
                            filtering_category = filtering_category,
                            filtering_value = filtering_value,
                            alpha_sign = as.numeric(alpha_sign),
                            rho_minimum = as.numeric(rho_minimum),
                            output_file_name_net = output_file_name_net,
                            network_categories_output = paste0(output_path, "/", community, "_", timepoints[[k]], "_when_pairwise_is_", consortia[[i]], "_", timepoints[[k]], "_and_",consortia[[j]],"_", timepoints[[k]], "_network_pathway_IDs"),
                            ref_value = ref_value)
                        
                        taxa_to_pathway_inocula = mapping_taxa_to_pathway(picrust2_file_strat = picrust2_file_strat, asv_to_taxa = asv_to_taxa, consortia = community)
                        
                        uncollapsed_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/", community, "_inocula_output/", community, "_pathway_out_contrib/", community, "_picrust2_pathway_with_taxa_names_uncollapsed.tsv")
                        
                        #encoded_matrix = make_encoded_matrix(per_seq_contrib = collapsed_file, asv_to_taxa_df = asv_to_taxa, taxa_list = rownames(enriched_clean_correlation_test), 
                        #                    pathway_list = colnames(enriched_clean_correlation_test), min_copy = 1e-9)
                        
                        heatmap_plotting(input_data_name = output_file_name_heat, matrix_to_plot = enriched_clean_correlation_test, 
                                            per_sequence_contribution_file = uncollapsed_file)  
                        # Inoc_file is taxa as rows and samples as columns
                        # Pathway_inocula_merged_metagenome has columns as pathways and samples as rows)
                        # asv_to_taxa has ASVs an taxa in each column.
                        # PICRUSt2 
                       #make_network(pathway_file_name = pathway_file,
                       #     taxa_file_name = taxa_file,
                       #     maaslin_file_name = maaslin_filt_file,
                       #     metadata_file_name = metadata_file,
                       #     output_path = output_path,
                       #     filtering_category = filtering_category,
                        #    filtering_value = filtering_value,
                        #    rho_minimum = as.numeric(rho_minimum),
                        #    community = community,
                        #    timepoint = timepoints[[k]],
                        #    clean_correlation_matrix_file = enriched_clean_correlation_test,
                        #    ref_value = ref_value,
                        #    output_file = output_file_name_net)
                    }
                }
            }
        }
    }
    else{
        which_correlation = "pathway" #Options are "taxa" or "pathways"
        for(i in seq(consortia)){
            for(j in seq(timepoints)){
                print("-----------------------------------")
                print(paste0("PROCESSING CONSORTIA ", consortia[[i]], " AT TIMEPOINT ", timepoints[[j]], " FOR SELF-COMPARISON"))
                taxa_file = paste0("/Users/danielcm/Desktop/diammatics/T1D/Phyloseq/ps_",tolower(consortia[[i]]),"_final.csv")
                self_clean_correlation_test = self_normalizing_and_filtering(
                    pathway_file = pathway_file,
                    taxa_file = taxa_file,
                    metadata_file = metadata_file,
                    output_path = paste0(heatmap_output_path, "/", consortia[[i]]),
                    filtering_category = filtering_category,
                    filtering_value = paste0(consortia[[i]], "_", timepoints[[j]]),
                    alpha_sign = as.numeric(alpha_sign),
                    rho_minimum = as.numeric(rho_minimum),
                    which_correlation = which_correlation)
        output_file_name_heat = paste0(heatmap_output_path, consortia[[i]],"_", timepoints[[j]], "_self_",which_correlation,"_correlation_analysis")
        uncollapsed_file_self = paste0("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/", consortia[[i]], "_inocula_output/", consortia[[i]], "_pathway_out_contrib/", consortia[[i]], "_picrust2_pathway_with_taxa_names_uncollapsed.tsv")
        heatmap_plotting(input_data_name = output_file_name_heat, matrix_to_plot = self_clean_correlation_test, per_sequence_contribution_file = uncollapsed_file_self)
            }
        }
    }
}

global()

library(tidyr)


