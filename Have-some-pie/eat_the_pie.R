# Daniel Castaneda Mogollon, PhD
# October 2nd, 2025
# Purpose: This script shows a workflow where it merges the taxonomy and the PICRUSt2 features to
# create a single plot showing the relationship amongs them.

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(compositions)

taxa_matrix = read.table("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/ps_weeks_final.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
pathway_matrix = read.table("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/Pathway_merged_metagenome.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata = read.table("/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)

subset_and_transpose = function(filter_to_use, subset, metadata) {
    samples_to_keep = rownames(metadata %>% filter(!!sym(filter_to_use) %in% subset))
    return(samples_to_keep)
}

samples_to_keep = subset_and_transpose("Week_and_consortia", c("S2_w5w6", "S2_w9w10"), metadata)
taxa_subset = taxa_matrix[,samples_to_keep]
pathway_subset = pathway_matrix[, samples_to_keep]
dim(taxa_subset)
dim(pathway_subset)


taxa_subset = t(taxa_subset)
pathway_subset = t(pathway_subset)

# Hellinger normalization for both matrices (after transposing)
taxa_subset_hel <- decostand(taxa_subset, method = "hellinger")
pathway_subset_hel <- decostand(pathway_subset, method = "hellinger")

# Calculate Bray-Curtis distances on Hellinger-normalized data
taxa_dist = vegdist(taxa_subset_hel, method = "bray")
pathway_dist = vegdist(pathway_subset_hel, method = "bray")

# PCoA
taxa_pcoa = cmdscale(taxa_dist, k = 2, eig = TRUE)
pathway_pcoa = cmdscale(pathway_dist, k = 2, eig = TRUE)

proc = procrustes(taxa_pcoa$points, pathway_pcoa$points)
proc_test = protest(taxa_pcoa$points, pathway_pcoa$points)
proc_test

taxa_scores = as.data.frame(taxa_pcoa$points)
taxa_scores$ID = rownames(taxa_scores)
taxa_scores = left_join(taxa_scores, metadata, by="ID")

pathway_scores <- as.data.frame(pathway_pcoa$points)
pathway_scores$ID <- rownames(pathway_scores)
pathway_scores <- left_join(pathway_scores, metadata, by="ID")

colnames(taxa_scores)[1:2] <- c("Dim1", "Dim2")
colnames(pathway_scores)[1:2] <- c("Dim1", "Dim2")

# Extract feature loadings (species scores) for taxa and pathways
taxa_fit = envfit(taxa_pcoa$points, taxa_subset_hel, permutations = 999)
taxa_vectors = as.data.frame(scores(taxa_fit, display="vectors"))
taxa_vectors$feature = rownames(taxa_vectors)
colnames(taxa_vectors)[1:2] <- c("Dim1", "Dim2")
taxa_vectors$vec_len = sqrt(taxa_vectors$Dim1^2 + taxa_vectors$Dim2^2)
top_taxa = taxa_vectors %>% arrange(desc(vec_len)) %>% head(7)

pathway_fit = envfit(pathway_pcoa$points, pathway_subset_hel, permutations = 999)
pathway_vectors = as.data.frame(scores(pathway_fit, display="vectors"))
pathway_vectors$feature = rownames(pathway_vectors)
colnames(pathway_vectors)[1:2] <- c("Dim1", "Dim2")
pathway_vectors$vec_len = sqrt(pathway_vectors$Dim1^2 + pathway_vectors$Dim2^2)
top_pathways = pathway_vectors %>% arrange(desc(vec_len)) %>% head(7)

# Add a Group column
taxa_scores$Group <- "Taxa"
pathway_scores$Group <- "Pathways"

ggplot() + 
  geom_point(
    data=taxa_scores,
    aes(Dim1, Dim2, fill=Merged_weeks, shape=Group),
    color="black", size=8, show.legend = TRUE
  ) +
  geom_point(
    data=pathway_scores,
    aes(Dim1, Dim2, fill=Merged_weeks, shape=Group),
    color="black", size=8, show.legend = TRUE
  ) +
  geom_segment(data=top_taxa, aes(x=0, y=0, xend=Dim1, yend=Dim2), 
               arrow=arrow(length=unit(0.3,"cm")), color="blue") +
  geom_label_repel(data=top_taxa, aes(x=Dim1, y=Dim2, label=feature), color="black", fill="white", size=8) +
  geom_segment(data=top_pathways, aes(x=0, y=0, xend=Dim1, yend=Dim2), 
               arrow=arrow(length=unit(0.3,"cm")), color="red") +
  geom_label_repel(data=top_pathways, aes(x=Dim1, y=Dim2, label=feature), color="black", fill="white", size=8) +
  scale_shape_manual(
    values = c("Taxa"=21,"Pathways" = 22)
  ) +
  scale_fill_manual(
    values = c(
      "w5w6" = "orange",
      "w9w10"  = "blue",
      "S5"  = "#7570b3",
      "NS6" = "#e7298a"
    )
  ) +
  guides(
    fill = guide_legend(title = "Merged_weeks", override.aes = list(shape = 21, size = 8, color = "black")),
    shape = guide_legend(title = "Group")
  ) +
  theme(
    text = element_text(size=25),
    axis.title = element_text(size=45),
    legend.text = element_text(size=25),
    axis.text = element_text(size = 35)
  )

# Optionally, remove or comment out base R biplots (they use PCA, not PCoA)
# biplot(taxa_pca, scaling=2)
# biplot(pathway_pca, scaling=2)
