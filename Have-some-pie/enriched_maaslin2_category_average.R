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


assigning_categories_to_enriched_pairwise_comparisons = function(cons1, cons2, comparison, communities_comp,subanalysis){
    if(comparison == "community"){
      df_enriched = read.csv(paste0(maaslin2_path,comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_","w5","_vs_",communities_comp[2],"_","w9w10","_ref_Week_and_Consortium,",communities_comp[2],"_","w9w10","/all_results.tsv"), sep = "\t")
      f_out = paste0(maaslin2_path,comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_","w5","_vs_",communities_comp[2],"_","w9w10","_ref_Week_and_Consortium,",communities_comp[2],"_","w9w10","/",communities_comp[1],"_","w5","_vs_",communities_comp[2],"_","w9w10","_enriched_pathways_with_categories.csv")

    }
    else{
        df_enriched = read.csv(paste0(maaslin2_path,comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/all_results.tsv"), sep = "\t")
        f_out = paste0(maaslin2_path,comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_enriched_pathways_with_categories.csv")
    }
    df_categories = read.csv("MetaCyc/Master_Files/Master_Metacyc_pathway_file_with_categories.tsv", sep = "\t")

    if(!("NegLog10_qval" %in% colnames(df_enriched))){
    print("NegLog10_qval column is not there, adding it now...")
    df_enriched$NegLog10_qval = -log10(df_enriched$qval)
    }

    significant_qval = 0.05
    significant_enriched = 1.0

    # Filter for significant pathways
    df_enriched_pwys = df_enriched[df_enriched$qval < significant_qval & abs(df_enriched$coef) > significant_enriched,] 
    df_enriched_pwys$feature = gsub(".", "-", df_enriched_pwys$feature, fixed = TRUE)

    level1_list = c("Level_1","Level_1.1","Level_1.2","Level_1.3") # Removing level 1 and level
    level2_list = c("Level_2","Level_2.1","Level_2.2","Level_2.3")
    
    all_levels = c(level1_list, level2_list)

    for(level in all_levels){
        print(level)
        df_enriched_pwys[[level]] = df_categories[[level]][match(df_enriched_pwys$feature, df_categories$PathwayID)] # match the pathway names in the enriched df to the master categories df and get the category for each pathway, for each level
        print(df_enriched_pwys[[level]])
    }

    df_enriched_pwys$pwy_name = df_categories$Names[match(df_enriched_pwys$feature, df_categories$PathwayID)]
    df_enriched_pwys$group = sapply(df_enriched_pwys$coef, function(x) ifelse(x > 0, paste0("Enriched in ", cons1), paste0("Enriched in ", cons2)))
    df_enriched_pwys$summary = df_categories$Summary[match(df_enriched_pwys$feature, df_categories$PathwayID)]

    print(paste0("Number of enriched pathways before filtering for unassigned categories: ", nrow(df_enriched_pwys)))
    #Remove # in command below if including level 1 pathways
    #df_enriched_pwys = df_enriched_pwys[!is.na(df_enriched_pwys$Level_1),] # filter out pathways that don't have an assigned category at level 1, since those are the ones we want to summarize in the boxplots
    print(paste0("Number of enriched pathways with assigned categories: ", nrow(df_enriched_pwys)))

    
    print(head(df_enriched_pwys[,1:16]))
    write.csv(df_enriched_pwys, f_out, row.names = FALSE)


    return(df_enriched_pwys)
}


making_figure = function(cons1, cons2, comparison, subanalysis){
  communities_comp = c(cons1, cons2)
  df_enriched_data = assigning_categories_to_enriched_pairwise_comparisons(cons1, cons2, comparison, communities_comp,subanalysis)
  
  unique_level1 = c(unique(df_enriched_data$Level_1), unique(df_enriched_data$Level_1.1), unique(df_enriched_data$Level_1.2), unique(df_enriched_data$Level_1.3))
  unique_level1 = unique(unique_level1[!is.na(unique_level1)])
  level1_list = c("Level_1","Level_1.1","Level_1.2","Level_1.3")
  level2_list = c("Level_2","Level_2.1","Level_2.2","Level_2.3")

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
  #for(category in level1_categories){
  #    seen_rows = c()
  #    values = c()
  #    for(levels in level1_list){
  #        matched_rows = which(!is.na(df_enriched_data[[levels]]) & df_enriched_data[[levels]] == category)
  #        new_rows = setdiff(matched_rows, seen_rows)  # avoid counting same pathway twice
  #        values = c(values, df_enriched_data$coef[new_rows])
  #        seen_rows = c(seen_rows, new_rows)
  #    }
  #    fold_values[[category]] = values
  #}

  #names(sort(fold_values[[category]], decreasing = TRUE))
  names(sort(fold_values2[[category]], decreasing = TRUE))

  #fold_values_all_levels = c(fold_values, fold_values2)
  fold_values_all_levels = fold_values2
  #print(max.levels=120,fold_values_all_levels)

  if(comparison == "community"){
    fold_tbl = enframe(fold_values_all_levels, name = "Category", value = "coef") %>%
    unnest(coef) %>%
    mutate(group = ifelse(coef > 0, paste0(communities_comp[1], " w5"), paste0(communities_comp[2], " w9w10")))
    fig_file = (paste0("/Users/danielcm/Desktop/SickKids/Maaslin2.6/",comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_","w5","_vs_",communities_comp[2],"_","w9w10","_ref_Week_and_Consortium,",communities_comp[2],"_","w9w10","/",communities_comp[1],"_","w5","_vs_",communities_comp[2],"_","w9w10","_all_enriched_category2.png"))
    group_levels = c(
      paste0(communities_comp[2], " w9w10"),
      paste0(communities_comp[1], " w5")
    )
    }

  else{
    fold_tbl <- enframe(fold_values_all_levels, name = "Category", value = "coef") %>%
    unnest(coef) %>%
    mutate(group = ifelse(coef > 0, paste0(communities_comp[1], " ", comparison), paste0(communities_comp[2], " ", comparison)))
    fig_file = (paste0("/Users/danielcm/Desktop/SickKids/Maaslin2.6/",comparison,"/PWY/",subanalysis,"/maaslin2_Week_and_Consortium_pathway_",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_ref_Week_and_Consortium,",communities_comp[2],"_",comparison,"/",communities_comp[1],"_",comparison,"_vs_",communities_comp[2],"_",comparison,"_all_enriched_category2.png"))
    group_levels = c(
      paste0(communities_comp[2], " ", comparison),
      paste0(communities_comp[1], " ", comparison)
    )
    }
  
    category_order = c(level1_categories, setdiff(level2_categories, level1_categories))
    fold_tbl$Category = factor(fold_tbl$Category, levels = rev(category_order))
    fold_tbl$group = factor(fold_tbl$group, levels = group_levels)

    boxplot_tbl = fold_tbl %>%
      group_by(Category, group) %>%
      filter(n() > 3) %>%
      ungroup()

    print(sort(fold_tbl[[category]], decreasing = TRUE), max.levels = 120)

    p = ggplot(fold_tbl, aes(x = coef, y = Category, fill = group)) +
    geom_boxplot(
      data = boxplot_tbl,
      aes(fill = group),
      outlier.shape = NA,
      alpha = 0.30,
      position = position_dodge(width = 0.75),
      outlier.stroke = 0.8
      ) +
      geom_point(
      aes(fill = group),
      shape = 21,
      size = 3,
      color = "black",
      stroke = 1,
      position = position_dodge(width = 0.75)
    ) +

    guides(fill = guide_legend(reverse = FALSE)) +

    scale_fill_manual(
      values = setNames(
        c("#F4B9C1","#36753B"),
        group_levels
      ),
      breaks = group_levels
    ) +
    theme_bw() +
    labs(
      x = "Log2(FC)",
    ) +
    theme(axis.text.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 20)) +

      scale_x_continuous(limits = c(-13, 13))

      ggsave(
        filename = fig_file,
        plot = p,
        width = 11,
        height = 6,
        #height = 11, (leave this 11 for cat1 and cat2)
        dpi = 600
      )
  print(p)
  return()
}
  # Summarize explicitly (avoids stat_summary quirks)
subanalysis_cat = c("community_contribution","microbe_contribution")
for(subanalysis in subanalysis_cat){
  making_figure("NS1", "S2", "w9w10", subanalysis)
  making_figure("NS6", "S5", "w9w10", subanalysis)
  making_figure("NS1", "S2", "w5", subanalysis)
  making_figure("NS6", "S5", "w5", subanalysis)
  making_figure("NS1", "NS1", "community", subanalysis)
  df_enriched_pwys = assigning_categories_to_enriched_pairwise_comparisons("NS1", "S2", "w9w10", c("NS1", "S2"), subanalysis)
}
