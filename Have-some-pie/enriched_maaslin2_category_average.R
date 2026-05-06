# Daniel Castaneda Mogollon, PhD
# February 22th, 2026
# This script, unlike 'Have-some-pie', will get the average per category level of all
# enriched pathways when doing a pairwise comparison with Maaslin2. 

# This script takes two inputs: the Maaslin2 all_results.tsv from the NS1 w9w10 vs S2 w9w10 pairwise pathway comparison, and a master MetaCyc reference file that maps pathway IDs to hierarchical categories (Levels 1–2) 
# and human-readable names. It filters the Maaslin2 results to pathways with qval < 0.05 and |coef| > 1, annotates each passing pathway with 
# its MetaCyc category hierarchy and enrichment group (NS1 or S2, based on the sign of the coefficient), and writes the annotated table to a CSV.
# It then groups the Maaslin2 coefficients by Level 2 MetaCyc category and produces a horizontal boxplot with individual points overlaid, saved as a 600 DPI PNG, 
# showing the distribution of fold-changes per category colored by which consortium was enriched.

# Load libraries ----------------------------------------------------------
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
# initializing()

master_path = "/Users/danielcm/Desktop/Sickkids/"
setwd(master_path)
maaslin2_path = paste0(master_path, "Maaslin2.6/")


assigning_categories_to_enriched_pairwise_comparisons = function(cons1, cons2, comparison, communities_comp){

    df_categories = read.csv("MetaCyc/Master_Files/Master_Metacyc_pathway_file_with_categories.tsv", sep = "\t")
    df_enriched = read.csv(paste0(maaslin2_path,comparison,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/all_results.tsv"), sep = "\t")

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
        df_enriched_pwys[[level]] = df_categories[[level]][match(df_enriched_pwys$feature, df_categories$Pathways)] # match the pathway names in the enriched df to the master categories df and get the category for each pathway, for each level
    }

    df_enriched_pwys$pwy_name = df_categories$Names[match(df_enriched_pwys$feature, df_categories$Pathways)]
    df_enriched_pwys$group = sapply(df_enriched_pwys$coef, function(x) ifelse(x > 0, paste0("Enriched in ", cons1), paste0("Enriched in ", cons2)))
    df_enriched_pwys$summary = df_categories$Summary[match(df_enriched_pwys$feature, df_categories$Pathways)]

    print(paste0("Number of enriched pathways before filtering for unassigned categories: ", nrow(df_enriched_pwys)))
    df_enriched_pwys = df_enriched_pwys[!is.na(df_enriched_pwys$Level.1),] # filter out pathways that don't have an assigned category at level 1, since those are the ones we want to summarize in the boxplots
    print(paste0("Number of enriched pathways with assigned categories: ", nrow(df_enriched_pwys)))

    write.csv(df_enriched_pwys, paste0(maaslin2_path,comparison,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_enriched_pathways_with_categories.csv"), row.names = FALSE)
    #View(df_enriched_pwys)

    return(df_enriched_pwys)
}




making_figure = function(cons1, cons2, comparison){
  communities_comp = c(cons1, cons2)
  df_enriched_data = assigning_categories_to_enriched_pairwise_comparisons(cons1, cons2, comparison, communities_comp)

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
      seen_rows = c()
      values = c()
      for(levels in level2_list){
          matched_rows = which(!is.na(df_enriched_data[[levels]]) & df_enriched_data[[levels]] == category)
          new_rows = setdiff(matched_rows, seen_rows)  # avoid counting same pathway twice
          values = c(values, df_enriched_data$coef[new_rows])
          seen_rows = c(seen_rows, new_rows)
      }
      fold_values2[[category]] = values
  }
  for(category in level1_categories){
      seen_rows = c()
      values = c()
      for(levels in level1_list){
          matched_rows = which(!is.na(df_enriched_data[[levels]]) & df_enriched_data[[levels]] == category)
          new_rows = setdiff(matched_rows, seen_rows)  # avoid counting same pathway twice
          values = c(values, df_enriched_data$coef[new_rows])
          seen_rows = c(seen_rows, new_rows)
      }
      fold_values[[category]] = values
  }

  names(fold_values)
  names(fold_values2)

  level1_categories


  fold_tbl <- enframe(fold_values2, name = "Category", value = "coef") %>%
    unnest(coef) %>%
    mutate(group = ifelse(coef > 0, paste0(communities_comp[1], " ", comparison), paste0(communities_comp[2], " ", comparison)))


    p = ggplot(fold_tbl, aes(x = coef, y = Category, fill = group)) +
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
      position = position_dodge(width = 0.75)
    ) +

    scale_fill_manual(
      values = setNames(
        c("darkorange1", "darkturquoise"),
        c(
          paste0(communities_comp[1], " ", comparison),
          paste0(communities_comp[2], " ", comparison)
        )
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

      ggsave(
        filename = paste0("/Users/danielcm/Desktop/SickKids/Maaslin2.6/",comparison,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_all_enriched.png"),
        plot = p,
        width = 10,
        height = 6,
        dpi = 600
      )
  return()
}
  # Summarize explicitly (avoids stat_summary quirks)

making_figure("NS1", "S2", "w9w10")
