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

new_filtering_parameter = function(community,rule_asv_count, rule_asv_length, rule_min_prevalence){
    ASV_count = rule_asv_count
    ASV_length = rule_asv_length
    prevalence_threshold = rule_min_prevalence

    keep = c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len")
    grepping = grepl(paste0(community, "(?!\\d)"),colnames(master_ps), perl = TRUE)
    keep = c(keep, colnames(master_ps)[grepping])
    subset_file = master_ps[, keep]
    #print(ncol(subset_file)) #135 columns for NS1, or 130 if removing metadata and keep only samples and inocula
    #View(subset_file)
    ignore_cols = c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len")
    subset_file$sum = rowSums(subset_file[, !colnames(subset_file) %in% ignore_cols])
    
    #Applying rules
    ignore_cols = c(ignore_cols, "sum")
    #Rule # 1: ASV count > 300
    df_1 = subset_file[subset_file$sum>ASV_count,]
    print(paste0("Number of ASVs after applying count filter: ", nrow(df_1)))

    #Rule # 2: ASV length > 250
    df_2 = df_1[df_1$asv_len > ASV_length,]
    print(paste0("Number of ASVs after applying length filter: ", nrow(df_2)))

    #Rule # 3: ASV present in at least 2 samples
    num_of_samples = ncol(df_2) - length(ignore_cols)
    df_2$nonzero_count = rowSums(df_2[, !colnames(df_2) %in% ignore_cols] != 0) #matches Laura's
    min_samples = ceiling(num_of_samples*prevalence_threshold)

    df_3 = df_2[df_2$nonzero_count >= (min_samples),]

    print(paste0("There are a total of ",sum(df_2$nonzero_count >= (min_samples)), " ASVs that are present in at least ", min_samples, " samples"))
    print(paste0("Number of ASVs after applying sample presence filter: ", nrow(df_3)))

    df_3$false_positive = ifelse(grepl(community, df_3$expected_communities), "No", "Yes")

    return(df_3)
    #Rule # 4: ASV relative abundance must be at least 10%
}
    
#apply_the_filter("NS1","inocula")
x = new_filtering_parameter("NS1", 300, 250, 0.10)

colnames(x)
sub_x = x[,colnames(x) %in% c("asv_id", "genus_final", "curated_species_femmicro","expected_communities", "asv_len", "sum", "nonzero_count")]

write.csv(x, "NS1_ASV_filtering_results.csv", row.names = FALSE)
