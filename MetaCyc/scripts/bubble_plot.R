# Daniel Castaneda Mogollon, PhD
# February 22nd, 2024
# This script creates a bubble plot to visualize the contribution of different taxa to the enriched pathways identified in the Maaslin2 analysis. The size of the bubbles represents the proportion of mice in which the taxa contribute to the pathway, while the color intensity represents the relative contribution of the taxa to that pathway.

library(ggplot2)

bubble = function(maaslin_f,df_contr,saving_file,direction){
    direction = "left"
    if(direction == "left"){
        enriched_pathways = df_maaslin$feature[df_maaslin$coef > 1 & df_maaslin$qval < 0.05]

    }
    elseif(direction == "right"){
        enriched_pathways = df_maaslin$feature[df_maaslin$coef < -1 & df_maaslin$qval < 0.05]
    }
    else{
        print("Direction not valid")
    }
    df_maaslin = read.csv(maaslin_f, sep = "\t")
    df_contribution = read.csv(df_contr, sep = "\t")

    taxa_present = unique(df_contribution$Taxa)
    taxa_present

    #View(bacteria_abundance_subset_df)

    df_contribution_enriched = df_contribution[df_contribution$pathway %in% enriched_pathways, ]
    #print(ncol(df_contribution_enriched))
    df_contribution_enriched = df_contribution_enriched[, !grepl("week5", colnames(df_contribution_enriched))]
    #print(ncol(df_contribution_enriched))

    ncol(df_contribution_enriched)

    #View(df_contribution_enriched)
    #Sanity check

    #print(nrow(df_contribution))
    #print(nrow(df_contribution_enriched))

    contribution_df = data.frame(Taxa = character(), Pathway = character(), Contribution = numeric(), Name = character())

    View(contribution_df)
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
            path_name = df_maaslin[df_maaslin$feature == pathway,'Names']
            num_of_mice = length(abund_per_mouse)
            if (num_of_mice == 0){
                #print(paste0("No abundance for ", bacteria, " in pathway ", pathway))
                next
            }
            prop_mice = num_of_mice/(ncol(df_contribution_enriched))
            bact_contribution = sum(abund_per_mouse)/num_of_mice
            df_tmp = data.frame(Taxa = bacteria, Pathway = pathway, Contribution = bact_contribution, Mice = num_of_mice, Prop_mice = prop_mice, Name = path_name)
            contribution_df = rbind(contribution_df, df_tmp)

    }
    #print(num_of_mice)
    #print(ncol(df_contribution_enriched))
    }


    for(pathway in unique(contribution_df$Pathway)){
        df_pathway = contribution_df[contribution_df$Pathway == pathway, ]
        df_pathway$Rel_Contribution = df_pathway$Contribution/sum(df_pathway$Contribution)
        contribution_df[contribution_df$Pathway == pathway, "Rel_Contribution"] = df_pathway$Rel_Contribution
    }
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
    scale_size_continuous(name = "Proportion of Mice", range = c(2, 12)) +
    scale_fill_gradient2(low = "#ffffff", mid = "darkorange", high = "darkred" ,midpoint = 0.47, name = "Relative Contribution") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90")
    ) + xlab("") + ylab("")

    print(bubble_plot)

    ggsave(saving_file, 
        plot = bubble_plot, 
        width = 14, height = 12, dpi = 600)

}

pairwise = list("w5","w9w10","consortium")
cons = list("NS1","NS6","S2","S5")
visualize = list("ns1","s2","ns6","s5")
path = "/Users/danielcm/Desktop/SickKids/Maaslin2.4/"
path_pie = "/Users/danielcm/Desktop/SickKids/PICRUSt2.3/"

bubble(maaslin_f = paste0(path,"w9w10/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_NS6_w9w10_ref_Week_and_consortia,NS6_w9w10/all_resultswith_categories.tsv"), 
                            df_contr = paste0(path_pie,"ns1_output/ns1_pathway_inference/path_abun_strat_with_taxa.tsv"),
                            saving_file = paste0(path,"w9w10/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_NS6_w9w10_ref_Week_and_consortia,NS6_w9w10/bubble_plot.png"))

 

