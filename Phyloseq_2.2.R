# Daniel Casataneda Mogollon, PhD
# February 22nd, 2026
# This script is a newer version of the phyloseq_t1d script I have been running

# This will be used for the filtering of the ASV tables in the inocula, 364 samples, and mice


###INITIALIZING PATH AND FILE LOCATION
path = "/Users/danielcm/Desktop/SickKids/Phyloseq2/"
setwd(path)                                                                     #Setting the path of where Im working

lauras_master_file = "/Users/danielcm/Desktop/SickKids/Phyloseq2/FemMicro_final_364_collapsed_all_JULY_TO_SHARE_20260306.csv"
master_ps = read.csv(lauras_master_file, sep = ",")

species_NS1_design = c("Alistipes finegoldii","Akkermansia muciniphila","Alistipes onderdonkii",
"Bifidobacterium longum","Bifidobacterium longum","Bifidobacterium pseudocatenulatum",
"Bacteroides stercoris","Bacteroides thetaiotaomicron","Bacteroides uniformis",
"Phocaeicola vulgatus","Collinsella sp902362275","Enterocloster bolteae",
"Enterocloster citroniae","Sarcina perfringens","Enterococcus_B durans",
"Eggerthella lenta","Escherichia coli","Eisenbergiella tayi","Flavonifractor plautii","Flavonifractor plautii",
"Hungatella effluvii","Intestinimonas_A celer","Parabacteroides distasonis","Parabacteroides merdae",
"Staphylococcus hominis","Streptococcus sp001556435", "Sutterella wadsworthensis")
species_NS1_design = unique(species_NS1_design)

species_NS6_design = c("Akkermansia massiliensis","Alistipes finegoldii","Alistipes onderdonkii",
"Bacteroides stercoris","Bacteroides stercoris","Bacteroides uniformis",
"Bacteroides uniformis","Phocaeicola vulgatus","Bacteroides xylanisolvens",
"Bifidobacterium adolescentis","Collinsella sp902362275","Bifidobacterium longum",
"Enterocloster bolteae","Enterocloster sp005845215","Sarcina perfringens",
"Eggerthella lenta","Eggerthella lenta","Eisenbergiella porci",
"Eisenbergiella porci","Eisenbergiella tayi","Enterococcus faecalis",
"Enterococcus faecalis","Escherichia coli","Flavonifractor plautii",
"Hungatella effluvii","Hungatella hathewayi","Parabacteroides distasonis",
"Parabacteroides merdae","Pseudoflavonifractor_A sp022772585","Staphylococcus hominis",
"Streptococcus sp001556435","Sutterella wadsworthensis")
species_NS6_design = unique(species_NS6_design)

species_S2_design = c("Akkermansia muciniphila","Bacillus subtilis","Phocaeicola dorei",
"Phocaeicola vulgatus","Bifidobacterium dentium","Bifidobacterium longum",
"Enterocloster bolteae","Alitiscatomonas sp900066535","Enterocloster clostridioformis",
"Otoolea symbiosa","Eggerthella lenta","Eggerthella lenta",
"Eisenbergiella tayi","Eisenbergiella tayi","Enterococcus_A avium",
"Enterococcus_A avium","Escherichia coli","Faecalibacterium longum",
"Flavonifractor plautii","Holdemania filiformis","Hungatella effluvii",
"Parabacteroides distasonis","Hungatella effluvii","Ruminococcus_B gnavus",
"Staphylococcus capitis","Staphylococcus hominis","Collinsella sp902472315", "Alistipes finegoldii",
"Roseburia_C amylophila")
species_S2_design = unique(species_S2_design)

species_S5_design = c("Collinsella sp902472315","Akkermansia muciniphila","Alistipes finegoldii",
"Alistipes finegoldii","Cytobacillus undefined","Phocaeicola vulgatus","Bifidobacterium longum",
"Bifidobacterium pseudocatenulatum","Enterocloster bolteae",
"Enterocloster citroniae","Enterocloster clostridioformis","Enterocloster alcoholdehydrogenati",
"Eggerthella lenta","Eisenbergiella porci","Eisenbergiella tayi",
"Enterococcus_A avium","Escherichia coli","Faecalibacterium longum",
"Flavonifractor plautii","Holdemania filiformis","Hungatella effluvii","CHH4-2 sp018378255",
"Parabacteroides distasonis","Parabacteroides distasonis","Agathobacter faecis","Ruminococcus_B gnavus",
"Staphylococcus epidermidis")
species_S5_design = unique(species_S5_design)

filter = function(community, data_type){
    identifiers = c("asv_id", "genus_final", "curated_species_femmicro")
    
    if(data_type == "mice"){
        keep = grepl(paste0(community, "(?!\\d)"),colnames(master_ps), perl = TRUE) | colnames(master_ps) %in% identifiers #The regex next to community avoids extracting other stuff for the S5 community, i.e. S58, S51, etc.
        subset_file = master_ps[, keep]
        keep2 = grepl("week5", colnames(subset_file)) | grepl("week9", colnames(subset_file)) | grepl("week10", colnames(subset_file)) | colnames (subset_file) %in% identifiers
        subset_file_f = subset_file[, keep2]
        print(paste0("Analyzing a total of ", ncol(subset_file_f), " columns after filtering for weeks"))
    return(subset_file_f)
    }

    else if(data_type == "inocula"){
        keep = grepl(community, colnames(master_ps)) | colnames(master_ps) %in% identifiers
        subset_file = master_ps[, keep]
        keep2 = grepl("ino", colnames(subset_file)) | colnames(subset_file) %in% identifiers
        subset_file = subset_file[, keep2]
        print(paste0("Analyzing a total of ", ncol(subset_file)-3, " columns for inocula:"))
        print(colnames(subset_file)[!colnames(subset_file) %in% identifiers])
        return(subset_file)
    }
    else{
        print("Data type not recognized, please input 'mice' or 'inocula'")
        return(NULL)
    }
    
}

ASV_count_hard_filter = function(df,rule){
    count_cols <- !colnames(df) %in% c("asv_id", "genus_final", "curated_species_femmicro")
    df$sum = rowSums(df[, count_cols])
    df_f = df[df$sum > rule,]
    taxa_kept = as.character(unique(paste0(df_f$genus_final," ",df_f$curated_species_femmicro)))
    print(paste0("Number of ASVs before filtering: ", nrow(df)))
    print(paste0("Number of ASVs after filtering: ", nrow(df_f)))

    return(taxa_kept)
}

true_positive_analysis = function(community,taxa_after_filter){
    if(community == "NS1"){
        design_taxa = species_NS1_design
    } else if (community == "NS6"){
        design_taxa = species_NS6_design
    } else if (community == "S2"){
        design_taxa = species_S2_design
    } else if (community == "S5"){
        design_taxa = species_S5_design
    } else {
        print("Community not recognized, please input NS1, NS6, S2, or S5")
        return(NULL)
    }

    TP = c()
    FP = c()
    FN = c()

    for(taxon in taxa_after_filter){ 
        if(taxon == "Enterocloster bolteae/citroniae"){ #This handles the Enterocloster discrepancy in FemMicro
            if("Enterocloster bolteae" %in% design_taxa) TP = c(TP, "Enterocloster bolteae")
            if("Enterocloster citroniae" %in% design_taxa) TP = c(TP, "Enterocloster citroniae")
        } 
        
        else if (taxon == "Enterocloster bolteae/clostridioformis"){
            if("Enterocloster bolteae" %in% design_taxa) TP = c(TP, "Enterocloster bolteae")
            if("Enterocloster clostridioformis" %in% design_taxa) TP = c(TP, "Enterocloster clostridioformis")
        } 
        
        else if (taxon == "Collinsella sp902362275/sp902472315"){
            if("Collinsella sp902362275" %in% design_taxa) TP = c(TP, "Collinsella sp902362275")
            if("Collinsella sp902472315" %in% design_taxa) TP = c(TP, "Collinsella sp902472315")
        }

        else if (taxon == "Parabacteroides distasonis/merdae"){
            if("Parabacteroides distasonis" %in% design_taxa) TP = c(TP, "Parabacteroides distasonis")
            if("Parabacteroides merdae" %in% design_taxa) TP = c(TP, "Parabacteroides merdae")
        } 
        
        else if (taxon == "Phocaeicola dorei/vulgatus"){
            if("Phocaeicola dorei" %in% design_taxa) TP = c(TP, "Phocaeicola dorei")
            if("Phocaeicola vulgatus" %in% design_taxa) TP = c(TP, "Phocaeicola vulgatus")
        }

        else{
            if(taxon %in% design_taxa){
                TP = c(TP, taxon)
            } else {
                FP = c(FP, taxon)
            }
        }
    }

    for(taxon in design_taxa){
        if(taxon %in% taxa_after_filter){
            next
        }
        else if(taxon == "Enterocloster bolteae" || taxon == "Enterocloster citroniae"){
            if("Enterocloster bolteae/citroniae" %in% taxa_after_filter | "Enterocloster bolteae/clostridioformis" %in% taxa_after_filter){
                 next
            } else {
                FN = c(FN, taxon)
            }
        } else if (taxon == "Enterocloster clostridioformis"){
            if("Enterocloster bolteae/clostridioformis" %in% taxa_after_filter){
                next
            } else {
                FN = c(FN, taxon)
            }
        } else if (taxon == "Collinsella sp902362275" || taxon == "Collinsella sp902472315"){
            if("Collinsella sp902362275/sp902472315" %in% taxa_after_filter){
                next
            } else {
                FN = c(FN, taxon)
            }
        } else if (taxon == "Parabacteroides distasonis" || taxon == "Parabacteroides merdae"){
            if("Parabacteroides distasonis/merdae" %in% taxa_after_filter){
                next
            } else {
                FN = c(FN, taxon)
            }
        } else if (taxon == "Phocaeicola dorei" || taxon == "Phocaeicola vulgatus"){
            if("Phocaeicola dorei/vulgatus" %in% taxa_after_filter){
                next
            } else {
                FN = c(FN, taxon)
            }
        }
         else {
             if(!(taxon %in% taxa_after_filter)){
                 FN = c(FN, taxon)
             }
        }
    }
    TP = unique(TP)
    FP = unique(FP)
    FN = unique(FN)

    print(paste0("True Positives: ", length(TP)))
    for(taxon in sort(TP)){
        cat(paste0(taxon,"\n"))
    }
    cat("\n")
    
    print(paste0("False Positives: ", length(FP)))
    for(taxon in sort(FP)){
        cat(paste0(taxon,"\n"))
    }
    cat("\n")

    print(paste0("False Negatives: ", length(FN)))
    for(taxon in sort(FN)){
        cat(paste0(taxon,"\n")) 
    }
    cat("\n")
    print(paste0("TP: ", length(TP), " FP: ", length(FP), " FN: ", length(FN)))

}

apply_the_filter = function(community,type){ #Community = "NS1", "NS6", "S2", or "S5". Type = "mice" or "inocula"
    df1_inocula = filter(community, type)
    taxa_kept_inocula = ASV_count_hard_filter(df1_inocula, 1000)
    analysis_inocula = true_positive_analysis(community, taxa_kept_inocula)
}


new_filtering_parameter = function(community, rule_asv_count, rule_asv_length, rule_min_prevalence){
    ASV_count = rule_asv_count
    ASV_length = rule_asv_length
    prevalence_threshold = rule_min_prevalence
    num_of_samples = 0

    keep = c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len")
    grepping = grepl(paste0(community, "(?!\\d)"),colnames(master_ps), perl = TRUE)
    keep = c(keep, colnames(master_ps)[grepping])
    subset_file = master_ps[, keep]
    if(community == "S5"){
        subset_file = subset(subset_file, select = -plate4_1186R_0_M_NS6_week5_S5_L0)
    }
    #print(ncol(subset_file)) #135 columns for NS1, or 130 if removing metadata and keep only samples and inocula
    #View(subset_file)
    ignore_cols = c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len")
    subset_file$sum = rowSums(subset_file[, !colnames(subset_file) %in% ignore_cols])
    
    invalid_species_vector = c()
    genus_and_species = paste0(subset_file$genus_final," ", subset_file$curated_species_femmicro)
    unique_genus_and_species = unique(genus_and_species)

    for(taxon in unique_genus_and_species){
        taxon_parts = strsplit(taxon, " ")[[1]]
        if(length(taxon_parts) < 2 || taxon_parts[length(taxon_parts)] == "NA"){
            invalid_species_vector = c(invalid_species_vector, taxon)
        }
    }

    some_taxa_with_no_species = head(unique(invalid_species_vector), 10)

    #print(paste0("Before filtering, there are a total of ", length(unique(invalid_species_vector)), " taxa with no species resolution"))
    #print(paste0("For example:" , paste(some_taxa_with_no_species, collapse = ", ")))

    #Applying rules
    ignore_cols = c(ignore_cols, "sum")
    #Rule # 1: ASV count > 300

    df_1 = subset_file[subset_file$sum>=ASV_count,]
    #print(paste0("Number of ASVs after applying count filter: ", nrow(df_1)))

    #Rule # 2: ASV length > 250
    df_2 = df_1[df_1$asv_len >= ASV_length,]
    #print(paste0("Number of ASVs after applying length filter: ", nrow(df_2)))

    #Rule # 3: ASV present in at least 2 samples
    num_of_samples = ncol(df_2) - length(ignore_cols)
    sample_cols = colnames(df_2)[!colnames(df_2) %in% ignore_cols]
    #print(sample_cols)
    df_2$nonzero_count = rowSums(df_2[, !colnames(df_2) %in% ignore_cols] != 0) #matches Laura's
    min_samples = ceiling(num_of_samples*prevalence_threshold)

    df_3 = df_2[df_2$nonzero_count >= (min_samples),]

    #print(paste0("There are a total of ",sum(df_2$nonzero_count >= (min_samples)), " ASVs that are present in at least ", min_samples, " samples"))
    #print(paste0("Number of ASVs after applying sample presence filter: ", nrow(df_3)))

    df_3$false_positive = ifelse(grepl(paste0(community, "(?!\\d)"), df_3$expected_communities, perl = TRUE), "No", "Yes")    
    sum_FP = sum(df_3$false_positive == "Yes")
    sum_TP = sum(df_3$false_positive == "No")
    PPV = (sum_TP/(sum_TP + sum_FP))*100
    asv_num = nrow(df_3)

    unique_taxa = unique(paste0(df_3$genus_final," ", df_3$curated_species_femmicro))
    num_unique_taxa = length(unique_taxa)

    genus_and_species_after_filtering = paste0(df_3$genus_final," ", df_3$curated_species_femmicro)
    invalid_species_vector_after_filtering = c()
    for(taxon in unique(genus_and_species_after_filtering)){
        taxon_parts = strsplit(taxon, " ")[[1]]
        if(length(taxon_parts) < 2 || taxon_parts[length(taxon_parts)]=="NA"){
            invalid_species_vector_after_filtering = c(invalid_species_vector_after_filtering, taxon)
        }
    }

    num_taxa_with_no_species_before_filter = length((invalid_species_vector))
    num_taxa_with_no_species = length(invalid_species_vector_after_filtering)

    #cat(paste0("Filter criteria: ASV count > ", ASV_count, ", ASV length > ", ASV_length, " ,ASV prevalence ", prevalence_threshold*100, 
    #"%, we get PPV:", PPV, "% TP:", sum_TP, " FP:", sum_FP, " Total ASVs: ", nrow(df_3), " Unique Taxa: ", num_unique_taxa, " Taxa with no species: ", num_taxa_with_no_species, "\n"))

    return(list(df_3, num_of_samples,min_samples,PPV, sum_TP, sum_FP, asv_num, num_unique_taxa, 
                num_taxa_with_no_species_before_filter,num_taxa_with_no_species))
    #Rule # 4: ASV relative abundance must be at least 10%
}
    
### This section onwards was done across the mouse data set (across all weeks and endpoints) with the inoculum samples as well (but no controls)
tunning_parameters = function(community){
    #This function will be used to tune the parameters of the filtering, i.e. the ASV count, ASV length, and prevalence threshold
    #It will return the PPV, TP, FP, and total ASVs for each combination of parameters
    df1 = data.frame()
    i = 1
    for(count in c(0,100,150,200,250,300,350,400,450,500,550,600)){ # 12 iterations
        for(min_len in c(0,200,250, 300, 350, 400)){ # 6 iterations
            for(prevalence in c(0,0.001,0.005,0.01,0.05, 0.10, 0.15, 0.20)){ # 8 iterations
                if(count == 0 && min_len == 0 && prevalence == 0 && i == 1){
                    tunning = new_filtering_parameter(community, 0, 0, 0) #Running just once with no filtering criteria, just to see how the raw processing looks like
                    i = i + 1
                    df1 = rbind(df1, c(0, 0, 0, unlist(tunning[2:10])))
                }
                else if(count == 0 | min_len == 0 | prevalence == 0){ #If the raw data has been calculated, then just jump to the next iteration set with new rules
                    next
                }
                else {
                    tunning = new_filtering_parameter(community, count, min_len, prevalence)
                    df1 = rbind(df1, c(count, min_len, prevalence, unlist(tunning[2:10]))) # Adding results as the iteration goes ...
                }
            }
        }
    }



    colnames(df1) = c("ASV_count", "ASV_length", "Prevalence_threshold", "Num_of_samples", "Min_samples_for_prevalence", paste0(community, "_PPV"), paste0(community, "_TP"), paste0(community, "_FP"), paste0(community, "_Total_ASVs"), 
                    paste0(community, "_Unique_taxa"), paste0(community, "_Taxa_with_no_species_before_filtering"), paste0(community, "_Taxa_with_no_species_after_filtering")) #Adding column names to the results data frame
                    
    write.csv(df1, paste0(community,"_filtering_tunning_results.csv"), row.names = FALSE) #Writing out the results across the community selected by the user
    return(df1)
}

merging_files_for_optimum_parameters = function(){ #Simply merges the files across all the iterations in each community, where the three parameters are acting as the key for merging the files
    file1 = read.csv("/Users/danielcm/Desktop/SickKids/Phyloseq2/NS1_filtering_tunning_results.csv")
    file2 = read.csv("/Users/danielcm/Desktop/SickKids/Phyloseq2/NS6_filtering_tunning_results.csv")
    file3 = read.csv("/Users/danielcm/Desktop/SickKids/Phyloseq2/S2_filtering_tunning_results.csv")
    file4 = read.csv("/Users/danielcm/Desktop/SickKids/Phyloseq2/S5_filtering_tunning_results.csv")

    file_f = merge(file1, file2, by = c("ASV_count", "ASV_length", "Prevalence_threshold"))
    file_f = merge(file_f, file3, by = c("ASV_count", "ASV_length", "Prevalence_threshold"))
    file_f = merge(file_f, file4, by = c("ASV_count", "ASV_length", "Prevalence_threshold"))
    write.csv(file_f, "NS1_S2_NS6_S5_filtering_tunning_results.csv", row.names = FALSE)
    return(file_f)
}

applying_best_filters_across_all = function(community,merged_file_tunning){
    ASV_count = 600
    ASV_length = 400
    prevalence_threshold = 0.20
    special_cols = c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len","asv_seq")
    grepping = grepl(paste0(community, "(?!\\d)"),colnames(master_ps), perl = TRUE)
    keep = c(special_cols, colnames(master_ps)[grepping])
    subset_file = master_ps[, keep]
    if(community == "S5"){
        subset_file = subset(subset_file, select = -plate4_1186R_0_M_NS6_week5_S5_L0)
    }

    num_of_samples = ncol(subset_file) - length(special_cols)
    sample_cols = colnames(subset_file)[!colnames(subset_file) %in% special_cols]

    subset_file$sum = rowSums(subset_file[, !colnames(subset_file) %in% special_cols])
    
    subset_file$nonzero_count = rowSums(subset_file[, !colnames(subset_file) %in% special_cols] != 0) #matches Laura's
    min_samples = ceiling(num_of_samples*prevalence_threshold)
    final_file = subset_file[subset_file$sum>=ASV_count,]
    final_file = final_file[final_file$asv_len >= ASV_length,]
    final_file = final_file[final_file$nonzero_count >= (min_samples),]
    write.csv(final_file, paste0(community,"_filtered_ASVs_count", ASV_count, "_len", ASV_length, "_prev", prevalence_threshold*100,".csv"), row.names = FALSE)
}

# Getting the results in each variable
results_NS1 = tunning_parameters("NS1")
results_S2 = tunning_parameters("S2")
results_S5 = tunning_parameters("S5")
results_NS6 = tunning_parameters("NS6")


merging_files_for_optimum_parameters()

# Applying the best filters selected by the user, the actual numbers are in the function, so they needed to be changed in there
applying_best_filters_across_all("NS1", merged_file_tunning)
applying_best_filters_across_all("NS6", merged_file_tunning)
applying_best_filters_across_all("S2", merged_file_tunning)
applying_best_filters_across_all("S5", merged_file_tunning)
