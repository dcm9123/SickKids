# Daniel Castaneda Mogollon, PhD
# February 22nd, 2024


df_maaslin = read.csv("/Users/danielcm/Desktop/SickKids/Maaslin2.3/pairwise_comparisons/w9w10/Pathway/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_consortia,S2_w9w10/NS1_w9w10_vs_S2_w9w10_enriched_pathways_with_categories.csv", sep = ",")
df_contribution = read.csv("/Users/danielcm/Desktop/SickKids/PICRUSt2/ns1_pathway_inference/path_abun_strat_with_taxa.tsv", sep = "\t")

enriched_pathways = df_maaslin$feature[df_maaslin$coef > 1 & df_maaslin$qval < 0.05]
taxa_present = unique(df_contribution$Taxa)

#View(bacteria_abundance_subset_df)

df_contribution_enriched = df_contribution[df_contribution$pathway %in% enriched_pathways, ]
print(ncol(df_contribution_enriched))
df_contribution_enriched = df_contribution_enriched[, !grepl("week5", colnames(df_contribution_enriched))]
print(ncol(df_contribution_enriched))

ncol(df_contribution_enriched)

#View(df_contribution_enriched)
#Sanity check

print(nrow(df_contribution))
print(nrow(df_contribution_enriched))

contribution_df = data.frame(Taxa = character(), Pathway = character(), Contribution = numeric())

for (bacteria in taxa_present) {
    for(pathway in enriched_pathways){
        df_bacteria = df_contribution_enriched[df_contribution_enriched$Taxa == bacteria, ]
        df_bacteria = df_bacteria[df_bacteria$pathway == pathway, ]
        if (nrow(df_bacteria) == 0){
            #print(paste0("No enriched pathways for ", bacteria))
            next
        } # bacteria has none of the enriched pathways
        abundance_matrix = as.matrix(df_bacteria[ , 3:(ncol(df_bacteria)-1)])
        abund_per_mouse = colSums(abundance_matrix)
        abund_per_mouse = abund_per_mouse[abund_per_mouse != 0]
        num_of_mice = length(abund_per_mouse)
        if (num_of_mice == 0){
            #print(paste0("No abundance for ", bacteria, " in pathway ", pathway))
            next
        }
        prop_mice = num_of_mice/(ncol(df_contribution_enriched))
        bact_contribution = sum(abund_per_mouse)/num_of_mice
        df_tmp = data.frame(Taxa = bacteria, Pathway = pathway, Contribution = bact_contribution, Mice = num_of_mice, Prop_mice = prop_mice)
        contribution_df = rbind(contribution_df, df_tmp)

        #print(paste0("Mean relative contribution of ", bacteria, " to ", pathway, ": ", bact_contribution))
    }
    print(num_of_mice)
    print(ncol(df_contribution_enriched))
}

View(contribution_df)

for(pathway in unique(contribution_df$Pathway)){
    df_pathway = contribution_df[contribution_df$Pathway == pathway, ]
    df_pathway$Rel_Contribution = df_pathway$Contribution/sum(df_pathway$Contribution)
    contribution_df[contribution_df$Pathway == pathway, "Rel_Contribution"] = df_pathway$Rel_Contribution
}



# Remove rows where species name is 'nan'
contribution_df = contribution_df[sapply(contribution_df$Taxa, function(x) {
    parts = strsplit(x, "_")[[1]]  # [[1]] extracts the vector from the list
    length(parts) < 2 || parts[2] != "nan"
}), ]
contribution_df$Taxa
View(contribution_df)

print(unique(contribution_df$Pathway))

contribution_df$Taxa_abbrev = sapply(contribution_df$Taxa, function(x) {
    parts = strsplit(x, "_")[[1]]  # [[1]] extracts the vector from the list
    if (length(parts) >= 2) {
        paste0(substr(parts[1], 1, 1), ". ", parts[2])
    } else {
        x
    }
})

library(ggplot2)

bubble_plot = ggplot(contribution_df, aes(x = Taxa_abbrev, y = Pathway, size = Prop_mice, fill = Rel_Contribution)) +
    geom_point(shape = 21, alpha = 0.8) +
    scale_size_continuous(name = "Proportion of Mice", range = c(2, 12)) +
    scale_fill_gradient(low = "yellow", high = "orangered", name = "Relative Contribution") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90")
    )

print(bubble_plot)

ggsave("/Users/danielcm/Desktop/SickKids/MetaCyc/figures/bubble_plot.png", 
       plot = bubble_plot, 
       width = 14, height = 16, dpi = 600)

