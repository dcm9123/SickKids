# Daniel Castaneda Mogollon, PhD
# October 5th, 2025
# Purpose: To evaluate different beta-diversity methods

library(ggplot2)
library(vegan)
library(phyloseq)


plotting_beta = function(df, meta, metatype, comparison,distance_type,method_type,color_values){
    meta = data.frame(ID=meta$ID, Subcommunity=meta$Subcommunity, Timepoint=meta$Timepoint, 
        Week_and_consortia = meta$Week_and_consortia, Sex_and_consortia = meta$Sex_and_consortia, 
        Sex_and_Timepoint = meta$Sex_and_timepoint)
    rownames(meta) = meta$ID  # SET THE ROW NAMES!
    ps1 = phyloseq(otu_table(otu_table(df, taxa_are_rows = TRUE)), sample_data(meta))
    ps1 = prune_samples(sample_data(ps1)[[metatype]] %in% comparison, ps1)  # Filter samples based on comparison
    ps1 = prune_samples(sample_sums(ps1) > 0, ps1)  # Remove samples with zero counts
    distance_matrix = vegdist(otu_table(ps1), method = distance_type)
    print(distance_matrix)
    ordination1 = ordinate(ps1, method = method_type, distance = distance_type)
    permanova = vegan::adonis2(distance_matrix ~ Subcommunity, data = sample_data(ps1))
    
    
    p1 = plot_ordination(ps1, ordination = ordination1, color = metatype) +
    geom_point(size = 4) +
    stat_ellipse(aes(group = .data[[metatype]]), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = color_values) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
  print(p1)
}


plotting_beta(df = read.table("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/Pathway_merged_metagenome.tsv", header = TRUE, sep = "\t", row.names=1),
               meta = read.table("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", 
                                 header = TRUE, sep = ",", row.names=1),
               metatype = "Subcommunity",
               comparison = c("NS1", "NS6", "S2", "S5"),
               distance_type = "bray",
               method_type = "NMDS",
               color_values = c("NS1" = "#E31A1C", "NS6" = "#1F78B4", 
                                "S2" = "#33A02C", "S5" = "#FF7F00"))




df = read.table("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/Pathway_merged_metagenome.tsv", header = TRUE, sep = "\t", row.names=1)
meta = read.table("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", 
                                 header = TRUE, sep = ",", row.names=1)
meta = data.frame(ID=meta$ID, Subcommunity=meta$Subcommunity, Timepoint=meta$Timepoint, 
        Week_and_consortia = meta$Week_and_consortia, Sex_and_consortia = meta$Sex_and_consortia, 
        Sex_and_Timepoint = meta$Sex_and_timepoint)
metatype = "Subcommunity"

    rownames(meta) = meta$ID  # SET THE ROW NAMES!
    ps1 = phyloseq(otu_table(otu_table(df, taxa_are_rows = TRUE)), sample_data(meta))
    ps1 = prune_samples(sample_data(ps1)[[metatype]] %in% comparison, ps1)  # Filter samples based on comparison
    ps1 = prune_samples(sample_sums(ps1) > 0, ps1)  # Remove samples with zero counts
    ps1
