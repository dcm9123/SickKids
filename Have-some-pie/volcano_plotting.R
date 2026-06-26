# Daniel Castaneda Mogollon, PhD
# November 6th, 2025
# Purpose: This script takes the output from Maaslin2 and the pairwise comparison from the data
# generated from FemMicro16S -> PICRUSt2 -> Maaslin2, and it creates a volcano plot, a vertical barplot
# showing the top 10 enriched categories on each side, and a donut chart showing the percentage
# of significant features per category.

#INPUT TYPES: dataframe from Maaslin2, dataframe with MetaCyc classification info, output path for the plots,
# category name for the left side, category name for the right side.

BiocManager::install("ggrepel")
BiocManager::install("ggpubr")

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggpubr)
library(beepr)


fetching_files = function(){
    df_maaslin2 = readline(prompt = "Enter the full path to the Maaslin2 output file. Must be a .tsv file e.g., /path/to/maaslin2_output/all_results.tsv: ")
    if(!file.exists(df_maaslin2)){
        stop("Maaslin2 output file does not exist. Please check the path and try again.")
    }
    pathway_classification = readline(prompt = "Enter the full path to the MetaCyc pathway classification file. Must be a .tsv file e.g., /path/to/maaslin2/metacyc_pathway_info.tsv: ")
    if(!file.exists(pathway_classification)){
        stop("MetaCyc pathway classification file does not exist. Please check the path and try again.")
    }
    
    category_left = readline(prompt = "Enter the specific name of the category on the left side of your plots (e.g., 'NS1 w9w10'): ")
    category_right = readline(prompt = "Enter the specific name of the category on the right side of your plots (e.g., 'S2 w9w10'): ")
    
    output_path = readline(prompt = "Enter the full path to the output directory e.g., /path/to/output_directory/: ")
    if(!dir.exists(output_path)){
        print("The output directory does not exist. Creating it now...")
        dir.create(output_path, recursive = TRUE)
    }
    return(list(df_maaslin2, pathway_classification, output_path, category_left, category_right))
}


enrichment_and_filtering = function(df_name, df_classification_name){

    df = read.table(df_name, header = TRUE, sep = "\t")                       #Reading the file provided by the user
    df_classification = read.csv(df_classification_name, header = FALSE, sep = "\t")          #Reading the classification file provided by the user

    df_classification = df_classification[,c(2,14,15)] #changed
    colnames(df_classification) = df_classification[1,]
    df_classification = df_classification[-1,] # Removing the first row
    df_classification[df_classification == ""] = NA

    print(dim(df_classification))
    print(dim(df))

    #print(head(df))

    df_merged = merge(df, df_classification, by.x = "feature", by.y = "PathwayID")          #Merging both dataframes to add classification info to the main dataframe and filling NAs
    print(head(df_merged))
    print(dim(df_merged))

    df = df[order(df$qval, decreasing = FALSE), ]                               #Ordering by q-value from zero to one
    df$neg_log10_qval = -log10(df$qval)                                         #Adding column of log10 qval conversion
    df$significant_right = ifelse(df$qval < 0.001 & df$coef > 1.0, "yes", "no") #Identifying significant points on the right side
    df$significant_left = ifelse(df$qval < 0.001 & df$coef < -1.0, "yes", "no") #Identifying significant points on the left side
    df$not_significant = ifelse(df$qval >= 0.001, "yes", "no")                  #Identifying non-significant points
    df$middle = ifelse(df$qval < 0.001 & abs(df$coef) <= 1.0, "yes", "no")      #Identifying points that are qval significant but with log2FC between -1 and 1
    total_significant = length(df$significant_left[df$significant_left == "yes"]) + 
                        length(df$significant_right[df$significant_right == "yes"]) #Counting total significant points from both sides

    return(list(df, total_significant)) #Returning the dataframe with new columns and the total significant points
}

volcano = function(df_volcano, category_left, category_right, output_path){
  axis_left = 0           #Creating variables to store the axis limits
  axis_right = 0          #Creating variables to store the axis limits

if(max(df_volcano$coef) > abs(min(df_volcano$coef))){ #Defining the range in a symmetrical way for the 'x' axis
  axis_right = round(max(df_volcano$coef)) + 1
  axis_left = (axis_right * -1)
} else {
  axis_left = round(min(df_volcano$coef)) - 1
  axis_right = abs(axis_left)
}

print(round(max(df_volcano$neg_log10_qval)))

x = ggplot(data = df_volcano, aes(x = coef, y = neg_log10_qval)) +      #Passing the df to ggplot
    geom_hline(yintercept = c(-log10(0.05), -log10(0.001)), linetype = "dashed", color = "black", linewidth = 0.8) + #Adding horizontal lines for qval thresholds
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.8) + #Adding vertical lines for log2FC thresholds

    #Adding jittered points to prevent overlapping
    geom_jitter(data = subset(df_volcano, not_significant == "yes"), shape = 21, color = "black", fill = "antiquewhite4", size = 4, width = 0.1, height = 0.05) +
    geom_jitter(data = subset(df_volcano, significant_right == "yes"), shape = 21, color = "black", fill = "#36753B", size = 4, width = 0.1, height = 0.05) +
    geom_jitter(data = subset(df_volcano, significant_left == "yes"), shape = 21, color = "black", fill = "#F4B9C1", size = 4, width = 0.1, height = 0.05) +
    geom_jitter(data = subset(df_volcano, middle == "yes"), shape = 21, color = "black", fill = "antiquewhite4", size = 4, width = 0.1, height = 0.05) +
    labs(
         x = paste0("← ", category_left, "      Log2(FC)        ", category_right, " →"), #Customizing x-axis label with categories
         y = "-Log10(q-val)") +
         theme(axis.text.x = element_text(size = 24), #Formatting text
               axis.text.y = element_text(size = 26),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 26),
               legend.title = element_text(size = 22),
               legend.text = element_text(size = 24),
               panel.background = element_rect(fill = "white", color = "black")) +
               scale_x_continuous(breaks = seq(axis_left, axis_right, by = 3), #Defining the range of the x axis
                                  limits = c(axis_left, axis_right)) +
               scale_y_continuous(limits = c(0, round(max(df_volcano$neg_log10_qval)) + 2),
                                  expand = c(0, 0)) # This removes padding below 0
    #geom_text_repel(data = subset(df_volcano, significant_right == "yes" & neg_log10_qval > 5), #Adding labels to significant points on the right side with high -log10 qval
    #                aes(label = feature))

    print(x)
    ggsave(paste0(output_path, "/volcano_plot.png"), plot = x, width = 7, height = 7, dpi = 600, bg = "white") #Saving the plot
}

vertical_barplot = function(df_volcano, output_path, category_left, category_right){
    df_volcano_sorted = df_volcano[order(df_volcano$coef, decreasing = TRUE), ] #Sorting the df by coef values
    rownames(df_volcano_sorted) = 1:nrow(df_volcano_sorted)                     #Adding a new index column
    top_ten_left = 0                #Variable for top left enriched categories
    top_ten_right = 0               #Variable for top right enriched categories
    
    top_ten_left = list(features = c(),coefs = c()) #Making a list for the name of the features and its log2FC
    top_ten_right = list(features = c(),coefs = c()) #Same but for top right side.
    for(i in (nrow(df_volcano_sorted)):(nrow(df_volcano_sorted)-9)){        #Loop that goes from the top 10 left enriched categories (most negative coef)
        if(df_volcano_sorted$significant_left[i] == "yes"){                 #If it's significant, then keep going
            top_ten_left[[1]] = c(top_ten_left[[1]], df_volcano_sorted$feature[i])
            top_ten_left[[2]] = c(top_ten_left[[2]], df_volcano_sorted$coef[i])
        }
        else{
            break  #If it's not significant, then break the loop
        }
    }

    #Does the same but for the right side
    for(j in 1:10){
        if(df_volcano_sorted$significant_right[j] == "yes"){
            top_ten_right[[1]] = c(top_ten_right[[1]], df_volcano_sorted$feature[j])
            top_ten_right[[2]] = c(top_ten_right[[2]], df_volcano_sorted$coef[j])
        }
        else{
            break
        }
    }

    order_index = order(top_ten_right$coefs, decreasing = FALSE) #Ordering the right side from lowest to highest log2FC for better visualization
    top_ten_right_sorted = list(features = top_ten_right$features[order_index], 
                            coefs = top_ten_right$coefs[order_index])
    #--SANITY CHECK --
    top_ten_right
    top_ten_right_sorted
    #--SANITY CHECK --

    df_left = data.frame(features = top_ten_left$features, coefs = top_ten_left$coefs, category = "left")
    df_right = data.frame(features = top_ten_right_sorted$features, coefs = top_ten_right_sorted$coefs, category = "right")

    #df_right
    #df_left
    df_both = rbind(df_left, df_right)
    write.table(x = df_both, file = paste0(output_path,"/vertical_plot_data.tsv"), sep = "\t") #Saving the dataframe used for plotting
    
    scale_left_side = round(min(df_both$coefs)) - 1
    scale_right_side = round(max(df_both$coefs)) + 1

    #df_both

    ggbarplot(df_both, x = "features", y = "coefs", color = "black") +
            coord_flip() + #Making the bars horizontal for better visualization
            scale_y_continuous(limits = c(scale_left_side, scale_right_side),
                            breaks = seq(scale_left_side, scale_right_side, by = 2)) +
            labs(y = paste0("← ", category_left, "  Log2(FC)   ", category_right, " →")) +
            geom_col(aes(fill = category), color = "black", width = 0.8) +
            scale_fill_manual(values = c("left" = "#F4B9C1", "right" = "#36753B")) +
            theme(axis.text.x = element_text(size = 28, face = "bold"),
                axis.text.y = element_text(size = 28, face = "bold"), 
                axis.title.x = element_text(size = 30, face = "bold"),
                axis.title.y = element_text(size = 0),
                legend.position = "none",
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 0),
                panel.background = element_rect(fill = "white", color = "black", size = 3)) +
                geom_hline(yintercept = 0, color = "black", linewidth = 1.5)
    
    ggsave(paste0(output_path,"/vertical_barplot.png"), width = 13, height = 9, dpi = 600, bg = "white")
    #Function to create vertical barplots showing the number of significant features per category
}

donut_chart = function(df_volcano, output_path, category_left, category_right){
    hole_size = 2
    group = c(category_left, category_right, "|Log2(FC)|<1 & q-val < 0.05", "Not Significant")
    values = c(length(df_volcano$significant_left[df_volcano$significant_left == "yes"]),
            length(df_volcano$significant_right[df_volcano$significant_right == "yes"]),
            length(df_volcano$middle[df_volcano$middle == "yes"]),
            length(df_volcano$not_significant[df_volcano$not_significant == "yes"]))

    total = length(df_volcano$feature)

    values_percentage = round((values / total) * 100, 2)
    sum_check = sum(values_percentage)
    if(sum_check!=100){
        max_value = which.max(values_percentage)
        values_percentage[max_value] = values_percentage[max_value]
        values_percentage[which.max(values_percentage)] = values_percentage[which.max(values_percentage)] - (sum_check - 100)
    }
    sum_check_m = sum(values_percentage)
    print(paste0("SANITY CHECK: THE ORIGINAL NUMBERS ADDED UP TO: ", sum_check, "% AND THE MODIFIED ONES ADDED UP TO: ", sum_check_m, "%"))

    donut_data = data.frame(category = paste(values_percentage,"% ",group), values = values, percentages = values_percentage, labels = values_percentage)

    colors2 = setNames(c("#F4B9C1",
                        "#36753B",
                        "white",
                        "antiquewhite4"), donut_data$category)

    names(colors2) = donut_data$category

    ggplot(donut_data, aes(x = hole_size, y = values, fill = category)) +
        geom_col(color = "black", linewidth = 1.5) + coord_polar(theta = "y") +
        scale_fill_manual(values = colors2) + 
        xlim(c(0.4, hole_size + 0.5)) +
        theme_void() +
        theme(legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 0),
            legend.key.size = unit(1, "cm"),
            plot.caption = element_text(size = 30, hjust = 0.5, face = "bold")) +
            labs(caption = paste("n = ", total))

    ggsave(paste0(output_path,"/donut_chart.png"), width = 9, height = 9, dpi = 600, bg = "white")
}

global = function(){
    path = "/Users/danielcm/Desktop/SickKids/Maaslin2.6/"
    group = c("communities","w5","w9w10","sex")
    for(g in group){
    files_to_process = list.files(path = paste0(path,g), pattern = "all_results.tsv", full.names = TRUE, recursive = TRUE)
        for(f in files_to_process){
            print(paste0("Processing file ",f))
            #"/Users/danielcm/Desktop/SickKids/Maaslin2.6/w9w10/PWY/microbe/maaslin2_Week_and_Consortium_pathway_NS1_w9w10_vs_NS6_w9w10_ref_Week_and_Consortium,NS6_w9w10/all_results.tsv"
            cat1 = strsplit(f,"/", fixed = TRUE)[[1]][10]
            print(cat1)
            cat1 = strsplit(cat1,"_",fixed = TRUE)[[1]][6:7]
            print(cat1)
            cat1 = paste(cat1[1],cat1[2],sep=' ')
            print(cat1)
            cat2 = strsplit(f,"/", fixed = TRUE)[[1]][10]
            cat2 = strsplit(cat2,"_",fixed = TRUE)[[1]][9:10]
            cat2 = paste(cat2[1],cat2[2],sep=" ")
            output_path = dirname(f)
            set_of_data = list(f,"/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/Master_Metacyc_pathway_file_with_categories.tsv",output_path,cat2,cat1)
            result = enrichment_and_filtering(df_name = set_of_data[[1]],
                                        df_classification = set_of_data[[2]])
            df_volcano_to_use = result[[1]]
            total_significant_points = result[[2]]

            volcano(df_volcano = df_volcano_to_use,
                    category_left = set_of_data[[4]],
                    category_right = set_of_data[[5]],
                    output_path = set_of_data[[3]])

            #vertical_barplot(df_volcano = df_volcano_to_use,
            #                output_path = set_of_data[[3]],
            #                category_left = set_of_data[[4]],
            #                category_right = set_of_data[[5]])

            donut_chart(df_volcano = df_volcano_to_use,
                        output_path = set_of_data[[3]],
                        category_left = set_of_data[[4]],
                        category_right = set_of_data[[5]])
            
            #beep(sound = 3)
        }
    }
}

global()




#Hello my name is Daniel and I am currently working on the volcano plot, vertical barplot,
# BUT I dontCheck()
#know shit
