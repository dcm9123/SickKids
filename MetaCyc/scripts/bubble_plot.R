# Daniel Castaneda Mogollon, PhD
# February 22nd, 2024
# This script creates a bubble plot to visualize the contribution of different taxa to the enriched pathways identified in the Maaslin2 analysis. The size of the bubbles represents the proportion of mice in which the taxa contribute to the pathway, while the color intensity represents the relative contribution of the taxa to that pathway.

library(ggplot2)

bubble = function(maaslin_f,df_contr,saving_file,direction){
    df_maaslin = read.csv(maaslin_f, sep = ",")
    df_contribution = read.csv(df_contr, sep = "\t")

    if(direction == "right"){
        enriched_pathways = df_maaslin$feature[df_maaslin$coef > 1 & df_maaslin$qval < 0.05]
    }

    else if(direction == "left"){
        enriched_pathways = df_maaslin$feature[df_maaslin$coef < -1 & df_maaslin$qval < 0.05]
    }

    else{
        stop("Direction not valid: use 'right' or 'left'")
    }


    taxa_present = unique(df_contribution$Taxa)
    #print(taxa_present)
    df_contribution_enriched = df_contribution[df_contribution$pathway %in% enriched_pathways, ]
    total_sample_cols = sum(!grepl("week5", colnames(df_contribution_enriched))) - 3  # total w9w10 sample columns, excluding metadata cols 1-3
    df_contribution_enriched = df_contribution_enriched[, !grepl("week5", colnames(df_contribution_enriched))]

    contribution_df = data.frame(Taxa = character(), Pathway = character(), Contribution = numeric(), Name = character())

    #View(contribution_df)
    for (bacteria in taxa_present) {
        for(pathway in enriched_pathways){
            df_bacteria = df_contribution_enriched[df_contribution_enriched$Taxa == bacteria, ]
            df_bacteria = df_bacteria[df_bacteria$pathway == pathway, ]
            if (nrow(df_bacteria) == 0){
                #print(paste0("No enriched pathways for ", bacteria))
                next
            } # bacteria has none of the enriched pathways
            abundance_matrix = as.matrix(df_bacteria[ , 4:(ncol(df_bacteria))])
            abund_per_mouse = colSums(abundance_matrix)
            abund_per_mouse = abund_per_mouse[abund_per_mouse != 0]
            path_name = df_maaslin[df_maaslin$feature == pathway,'pwy_name'] # get the human-readable pathway name from the maaslin df based on the pathway ID
            num_of_mice = length(abund_per_mouse)
            if (num_of_mice == 0){
                print(paste0("No abundance for ", bacteria, " in pathway ", pathway))
                next
            }
            prop_mice = num_of_mice/total_sample_cols

            bact_contribution = sum(abund_per_mouse)/num_of_mice

            # Troubleshooting
            #print(length(bacteria))
            #print(length(pathway))
            #print(length(bact_contribution))
            #print(length(num_of_mice))
            #print(length(prop_mice))
            #print(length(path_name))

            df_tmp = data.frame(Taxa = bacteria, Pathway = pathway, Contribution = bact_contribution, Mice = num_of_mice, Prop_mice = prop_mice, Name = path_name)
            #View(df_tmp)

            contribution_df = rbind(contribution_df, df_tmp)
        }
    }


    for(pathway in unique(contribution_df$Pathway)){
        df_pathway = contribution_df[contribution_df$Pathway == pathway, ]
        df_pathway$Rel_Contribution = df_pathway$Contribution/sum(df_pathway$Contribution)
        contribution_df[contribution_df$Pathway == pathway, "Rel_Contribution"] = df_pathway$Rel_Contribution
    }
    print(paste0("prop_mice range: ", round(min(contribution_df$Prop_mice),3), " - ", round(max(contribution_df$Prop_mice),3)))
    genus_and_species_name <- sapply(contribution_df$Taxa, function(x) {
        parts <- strsplit(x, "_")[[1]]
        if (length(parts) < 2) return(x)
            genus <- parts[1]
    # species: keep last entry if slash-separated
            species <- tail(strsplit(parts[2], "/")[[1]], 1)
        if (length(parts) >= 3) {
            paste(genus, species, parts[3])
        } else {
            paste(genus, species)
        }
    })

    contribution_df$Taxa_abbrev <- genus_and_species_name
    contribution_df <- contribution_df[
    !grepl("_nan$", contribution_df$Taxa),]


    bubble_plot = ggplot(contribution_df, aes(x = Taxa_abbrev, y = Pathway, size = Prop_mice, fill = Rel_Contribution)) +
    geom_point(shape = 21, alpha = 0.8) +
    scale_size_continuous(name = "Proportion of Mice", range = c(4, 18), limits = c(0.30, 1), breaks = c(0.30, 0.50, 0.70, 0.90, 1.00)) +
    scale_fill_gradient2(low = "#ffffff", mid = "darkorange", high = "darkred" ,midpoint = 0.47, name = "Relative Contribution",) +
    guides(fill = guide_colorbar(title.vjust = 2)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 26, face = "bold"),
        axis.text.y = element_text(size = 28, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28, face = "bold"),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(1.5, "cm"),
        panel.grid.major = element_line(color = "grey90")
    ) + xlab("") + ylab("")

    print(bubble_plot)

    ggsave(saving_file,
        plot = bubble_plot,
        width = 18, height = 14, dpi = 600)

}


bubble(maaslin_f = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_Consortium,S2_w9w10/NS1_w9w10_vs_S2_w9w10_enriched_pathways_with_categories.csv"),
                            df_contr = paste0(path_pie,"ns1_output/ns1_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_Consortium,S2_w9w10/ns1_w9w10_when_compared_s2_w9w10_bubble_plot.png"),
                            direction = "right") #NS1 enriched pathways and its contributing taxa in w9w10 comparison

bubble(maaslin_f = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_Consortium,S2_w9w10/NS1_w9w10_vs_S2_w9w10_enriched_pathways_with_categories.csv"),
                            df_contr = paste0(path_pie,"s2_output/s2_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_Consortium,S2_w9w10/s2_w9w10_when_compared_ns1_w9w10_bubble_plot.png"),
                            direction = "left") #S2 enriched pathways and its contributing taxa in w9w10 comparison

bubble(maaslin_f = paste0(path,"w5/maaslin2_Week_and_Consortium_pathway_NS1_w5w6_vs_S2_w5w6_ref_Week_and_Consortium,S2_w5w6/all_resultswith_categories.tsv"),
                            df_contr = paste0(path_pie,"ns1_output/ns1_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w5/maaslin2_Week_and_Consortium_pathway_NS1_w5w6_vs_S2_w5w6_ref_Week_and_Consortium,S2_w5w6/ns1_w5w6_when_compared_s2_w5w6_bubble_plot.png"),
                            direction = "right")

bubble(maaslin_f = paste0(path,"w5/maaslin2_Week_and_Consortium_pathway_NS1_w5w6_vs_S2_w5w6_ref_Week_and_Consortium,S2_w5w6/all_resultswith_categories.tsv"),
                            df_contr = paste0(path_pie,"s2_output/s2_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w5/maaslin2_Week_and_Consortium_pathway_NS1_w5w6_vs_S2_w5w6_ref_Week_and_Consortium,S2_w5w6/s2_w5w6_when_compared_ns1_w5w6_bubble_plot.png"),
                            direction = "left")

bubble(maaslin_f = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS6_w9w10_vs_S5_w9w10_ref_Week_and_Consortium,S5_w9w10/NS6_w9w10_vs_S5_w9w10_enriched_pathways_with_categories.csv"),
                            df_contr = paste0(path_pie,"ns6_output/ns6_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS6_w9w10_vs_S5_w9w10_ref_Week_and_Consortium,S5_w9w10/ns6_w9w10_when_compared_s5_w9w10_bubble_plot.png"),
                            direction = "right") #S5 enriched pathways and its contributing taxa in w9w10 comparison


bubble(maaslin_f = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS6_w9w10_vs_S5_w9w10_ref_Week_and_Consortium,S5_w9w10/NS6_w9w10_vs_S5_w9w10_enriched_pathways_with_categories.csv"),
                            df_contr = paste0(path_pie,"s5_output/s5_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w9w10/maaslin2_Week_and_Consortium_pathway_NS6_w9w10_vs_S5_w9w10_ref_Week_and_Consortium,S5_w9w10/s5_w9w10_when_compared_ns6_w9w10_bubble_plot.png"),
                            direction = "left") #S5 enriched pathways and its contributing taxa in w9w10 comparison


