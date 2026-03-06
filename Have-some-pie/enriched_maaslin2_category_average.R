# Daniel Castaneda Mogollon, PhD
# February 22th, 2026
# This script, unlike 'Have-some-pie' will get the average per category level of all
# enriched pathways when doing a pairwise comparison with Maaslin2. 

# Load libraries ----------------------------------------------------------
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
# initializing()

master_path = "/Users/danielcm/Desktop/Sickkids/"
maaslin2_path = paste0(master_path, "Maaslin2.4/")

assigning_categories_to_enriched_pairwise_comparisons = function(cons1, cons2){
    setwd(master_path)
    df_categories = read.csv("MetaCyc/Master_Files/Master_Metacyc_pathway_file_with_categories.tsv", sep = "\t")
    df_enriched = read.csv(paste0(maaslin2_path,"w9w10/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_consortia,S2_w9w10/all_results.tsv"), sep = "\t")

    if(!("NegLog10_qval" %in% colnames(df_enriched))){
    print("NegLog10_qval column is not there, adding it now...")
    df_enriched$NegLog10_qval = -log10(df_enriched$qval)
    }

    significant_qval = 0.05
    significant_enriched = 1.0

    # Filter for significant pathways
    df_enriched_pwys = df_enriched[df_enriched$qval < significant_qval & abs(df_enriched$coef) > significant_enriched,] 
    df_enriched_pwys$feature = gsub(".", "-", df_enriched_pwys$feature, fixed = TRUE)
    print(df_enriched_pwys$feature)


    level1_list = c("Level.1","Level.1.1","Level.1.2","Level.1.3")
    level2_list = c("Level.2","Level.2.1","Level.2.2","Level.2.3")
    all_levels = c(level1_list, level2_list)

    for(level in all_levels){
        df_enriched_pwys[[level]] = df_categories[[level]][match(df_enriched_pwys$feature, df_categories$Pathways)]
    }

    df_enriched_pwys$pwy_name = df_categories$Names[match(df_enriched_pwys$feature, df_categories$Pathways)]
    df_enriched_pwys$group = sapply(df_enriched_pwys$coef, function(x) ifelse(x > 0, paste0("Enriched in ", cons1), paste0("Enriched in ", cons2)))
    df_enriched_pwys$summary = df_categories$Summary[match(df_enriched_pwys$feature, df_categories$Pathways)]

    print(paste0("Number of enriched pathways before filtering for unassigned categories: ", nrow(df_enriched_pwys)))
    df_enriched_pwys = df_enriched_pwys[!is.na(df_enriched_pwys$Level.1),]
    print(paste0("Number of enriched pathways with assigned categories: ", nrow(df_enriched_pwys)))

    write.csv(df_enriched_pwys, paste0(maaslin2_path,"w9w10/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_consortia,S2_w9w10/NS1_w9w10_vs_S2_w9w10_enriched_pathways_with_categories.csv"), row.names = FALSE)
    #View(df_enriched_pwys)

    return(df_enriched_pwys)
}

plotting_data = function(df){

}

df_enriched_data = assigning_categories_to_enriched_pairwise_comparisons("NS1 w9w10", "S2 w9w10")

unique_level1 = c(unique(df_enriched_data$Level.1), unique(df_enriched_data$Level.1.1), unique(df_enriched_data$Level.1.2), unique(df_enriched_data$Level.1.3))
unique_level1 = unique(unique_level1[!is.na(unique_level1)])
level1_list = c("Level.1","Level.1.1","Level.1.2","Level.1.3")
level2_list = c("Level.2","Level.2.1","Level.2.2","Level.2.3")

level1_categories = c()
level2_categories = c()

level1_categories <- unique(unlist(df_enriched_data[level1_list]))
level1_categories <- level1_categories[!is.na(level1_categories) & level1_categories != ""]

level2_categories <- unique(unlist(df_enriched_data[level2_list]))
level2_categories <- level2_categories[!is.na(level2_categories) & level2_categories != ""]

fold_values = list()
fold_values2 = list()

for(category in level2_categories){
    values = c()
    for(levels in level2_list){
        values = c(values, df_enriched_data$coef[df_enriched_data[[levels]] == category])
    }
    fold_values2[[category]] = values
}
for(category in level1_categories){
    values = c()
    for(levels in level1_list){
        values = c(values, df_enriched_data$coef[df_enriched_data[[levels]] == category])
    }
    fold_values[[category]] = values
}

names(fold_values)
names(fold_values2)




fold_tbl <- enframe(fold_values2, name = "Category", value = "coef") %>%
  unnest(coef) %>%
  mutate(group = ifelse(coef > 1, "NS1 w9w10", "S2 w9w10"))


  ggplot(fold_tbl, aes(x = coef, y = Category, fill = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
    geom_point(
    aes(fill = group),
    shape = 21,
    size = 3,
    color = "black",
    stroke = 0.5,
    position = position_dodge(width = 0.5)
  ) +

  scale_fill_manual(
    values = c(
      "NS1 w9w10" = "darkorange1",
      "S2 w9w10" = "darkturquoise"
    )
  ) +
  theme_bw() +
  labs(
    x = "Log2(FC)",
  ) +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +

    scale_x_continuous(limits = c(-7.5, 7.5))
    
    ggsave("/Users/danielcm/Desktop/SickKids/Maaslin2.4/w9w10/maaslin2_Week_and_consortia_pathway_NS1_w9w10_vs_S2_w9w10_ref_Week_and_consortia,S2_w9w10/NS1_w9w10_vs_S2_w9w10_all_enriched.png", width = 10, height = 6, dpi = 600)

# Summarize explicitly (avoids stat_summary quirks)
