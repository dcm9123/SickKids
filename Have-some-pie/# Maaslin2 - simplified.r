# Maaslin2 - simplified
# April 2026

install.packages("beepr")
biocmanager::install("beepr")

library(Maaslin2)
library(beepr)

# Defining global variables
project_path = file.path(Sys.getenv("HOME"), "Desktop", "SickKids")
path_picrust2 = file.path(project_path, "PICRUSt2.6")
maaslin_path = file.path(project_path, "Maaslin2.6")
metadata_file = file.path(project_path, "Metadata", "Danska_diabetes_metadata364_20260409.csv")
normalization_method = "TSS"
transformation = "LOG"
cons1 = c("NS1","NS6","S2","S5")
cons2 = c("NS1","NS6","S2","S5")
metadata_id_to_match = "SRA_sample_name"

picrust_to_maaslin2<-function(file_to_use, variable, 
                              feature_to_use,random_input,reference_val,
                              subset,subset1,subset2, functions_to_analyze){
  if (!dir.exists(maaslin_path)) {
    dir.create(maaslin_path, recursive = TRUE, showWarnings = FALSE)
  }
  input_file = file.path(path_picrust2, file_to_use)
  if (!file.exists(input_file)) {
    stop("PICRUSt2 input file not found: ", input_file)
  }
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  setwd(path_picrust2)
  maaslin_input = read.table(input_file, header = TRUE, row.names = 1, sep='\t')
  metadata_input = read.csv(metadata_file, header = TRUE, row.names = metadata_id_to_match)
  fixed_effects_to_use = variable
  transformation_method = transformation
  normalization_to_use = normalization_method
  random_to_use = random_input
  reference_value = reference_val

  if(subset == TRUE){
    metadata = metadata_input[metadata_input[,variable] %in% c(subset1,subset2),]
    output_file = file.path(maaslin_path, paste0("maaslin2_",variable,"_",feature_to_use,"_",subset1,"_vs_",subset2,"_ref_",reference_value))
    Maaslin2(input_data = maaslin_input, input_metadata = metadata, output=output_file,
    fixed_effects = fixed_effects_to_use, transform = transformation_method,
    normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
    cores = 8, plot_scatter = FALSE)
  } 
  else {
    output_file = file.path(maaslin_path, paste0("maaslin2_",variable,"_",feature_to_use,"_",reference_value))
    Maaslin2(input_data = maaslin_input, input_metadata = metadata_input, output=output_file, # nolint
    fixed_effects = fixed_effects_to_use, transform = transformation_method, 
    normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
    cores = 6, plot_scatter = FALSE)
  }
  
  cat("The results were generated based on the following information: \n")
  cat(paste0("Input file: ",input_file,"\n"))
  cat(paste0("Metadata file: ",metadata_file,"\n"))
  cat(paste0("Fixed effect variable: ",variable,"\n"))
  cat(paste0("Random effect variable: ",random_input,"\n"))
  cat(paste0("Reference value: ",reference_val,"\n"))
  cat(paste0("Feature: ",feature_to_use,"\n"))
  cat(paste0("Normalization method: ",normalization_method,"\n"))
  cat(paste0("Transformation method: ",transformation,"\n"))
  cat(paste0("File written to: ",output_file,"\n"))
  cat("\n")
  cat("Finished!")
}

subcommunities_diff_time = function(functions_to_analyze, file_suffix_to_use){
  for (i in cons1[1:4]) # DONE
  {
    for (k in functions_to_analyze)
    {
      feature = k
      subs = paste0(i,"_w5_vs_",i,"_w9w10")
      picrust_to_maaslin2(variable = "Week_and_Consortium",
                          file_to_use = paste0(feature, file_suffix_to_use, ".tsv"),
                          feature_to_use = feature,
                          random_input = "Sex",
                          reference_val = paste0("Week_and_Consortium,",i,"_w9w10"),
                          subset = TRUE, subset1 = paste0(i,"_w5"), subset2 = paste0(i,"_w9w10"))
      beepr::beep(3)
      print(paste0("Finished ", subs, " ", feature))
    }
  }
}

time_w5_diff_subcommunities = function(functions_to_analyze, file_suffix_to_use){
### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W5### - DONE
func = functions_to_analyze
  for (i in cons1[1:4]){
    for (j in cons2[1:4]){
      if (j<=i){
        next()
      }
      else{
        for (k in func){
          feature = k
          subs = paste0(i,"_w5_vs_",j,"_w5")
          picrust_to_maaslin2( # nolint
                              variable = "Week_and_Consortium",
                              file_to_use = paste0(feature, file_suffix_to_use, ".tsv"),
                              feature_to_use = feature,
                              random_input = "Sex",
                              reference_val = paste0("Week_and_Consortium,",j,"_w5"),
                              subset = TRUE, subset1 = paste0(i,"_w5"), subset2 = paste0(j,"_w5"))
          beepr::beep(3)
          print(paste0("Finished ", subs, " ", feature))
        }
      }
    }
  }
}

time_w9w10_diff_subcommunities = function(functions_to_analyze, file_suffix_to_use){
  func = functions_to_analyze
  for (i in cons1[1:4]){
    for (j in cons2[1:4]){
      if (j<=i){
        next()
      }
      else{
        for (k in func){
          feature = k
          subs = paste0(i,"_w9w10_vs_",j,"_w9w10")
          picrust_to_maaslin2( # nolint
                              variable = "Week_and_Consortium",
                              file_to_use = paste0(feature, file_suffix_to_use, ".tsv"),
                              feature_to_use = feature,
                              random_input = "Sex",
                              reference_val = paste0("Week_and_Consortium,",j,"_w9w10"),
                              subset = TRUE, subset1 = paste0(i,"_w9w10"), subset2 = paste0(j,"_w9w10"))
          beepr::beep(3)
          print(paste0("Finished ", subs, " ", feature))
        }
      }
    }
  }
}

same_subcommunities_diff_sex = function(functions_to_analyze, file_suffix_to_use){
  for (subs in cons1[1:4]){
    for (k in functions_to_analyze){
      feature = k
      picrust_to_maaslin2( # nolint
                          variable = "Sex_and_Consortium",
                          file_to_use = paste0(feature, file_suffix_to_use, ".tsv"),
                          feature_to_use = feature,
                          random_input = NULL,
                          reference_val = paste0("Sex_and_Consortium,Male", "_", subs),
                          subset = TRUE, subset1 = paste0("Male", "_", subs), subset2 = paste0("Female", "_", subs))
      beepr::beep(3)
      print(paste0("Finished ", subs, " ", feature))
    }
  }
}

main_function = function(functions_to_analyze, file_suffix_to_use){
  subcommunities_diff_time(functions_to_analyze, file_suffix_to_use)
  time_w5_diff_subcommunities(functions_to_analyze, file_suffix_to_use)
  time_w9w10_diff_subcommunities(functions_to_analyze, file_suffix_to_use)
  same_subcommunities_diff_sex(functions_to_analyze, file_suffix_to_use)
}

fun_to_analyze = "pathway"
#fun_to_analyze = c("KO","EC","pathway")
main_function(functions_to_analyze = fun_to_analyze, file_suffix_to_use = "_merged_metagenome_strat_by_pwy_no_taxa")
#main_function(functions_to_analyze = fun_to_analyze, file_suffix_to_use = "_merged_metagenome")
beepr::beep(8)
