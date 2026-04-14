# Maaslin2 - simplified
# April 2026

library(Maaslin2)
library(beepr)

# Defining global variables
path_picrust2 = "/Users/danielcm/Desktop/Sickkids/PICRUSt2.6/"
maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.6/"
metadata_file = "/Users/danielcm/Desktop/SickKids/Metadata/Danska_diabetes_metadata364_20260409.csv"
normalization_method = "TSS"
transformation = "LOG"
cons1 = c("NS1","NS6","S2","S5")
cons2 = c("NS1","NS6","S2","S5")
metadata_id_to_match = "SRA_sample_name"

picrust_to_maaslin2<-function(file_to_use, variable, 
                              feature_to_use,random_input,reference_val,
                              subset,subset1,subset2){
  setwd(path_picrust2)
  maaslin_input = read.table(file_to_use, header = TRUE, row.names = 1, sep='\t')
  metadata_input = read.csv(metadata_file, header = TRUE, row.names = metadata_id_to_match)
  fixed_effects_to_use = variable
  transformation_method = transformation
  normalization_to_use = normalization_method
  random_to_use = random_input
  reference_value = reference_val
  if(subset == TRUE){
    metadata = metadata_input[metadata_input[,variable] %in% c(subset1,subset2),]
    output_file = paste0(maaslin_path,"maaslin2_",variable,"_",feature_to_use,"_",subset1,"_vs_",subset2,"_ref_",reference_value)
    Maaslin2(input_data = maaslin_input, input_metadata = metadata, output=output_file,
    fixed_effects = fixed_effects_to_use, transform = transformation_method,
    normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
    cores = 8, plot_scatter = FALSE)
  } else {
    output_file = paste0(maaslin_path,"maaslin2_",variable,"_",feature_to_use,"_",reference_value)
    Maaslin2(input_data = maaslin_input, input_metadata = metadata_input, output=output_file, # nolint
    fixed_effects = fixed_effects_to_use, transform = transformation_method, 
    normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
    cores = 6, plot_scatter = FALSE)
  }
  cat("The results were generated based on the following information: \n")
  cat(paste0("Input file: ",path_picrust2,file_to_use,"\n"))
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

subcommunities_diff_time = function(){
  for (i in cons1[1:4]) # DONE
  {
    for (k in c("EC","KO","pathway"))
    {
      feature = k
      subs = paste0(i,"_w5_vs_",i,"_w9w10")
      picrust_to_maaslin2(variable = "Week_and_Consortium",
                          file_to_use = paste0(feature,"_merged_metagenome.tsv"),
                          feature_to_use = feature,
                          random_input = "Sex",
                          reference_val = paste0("Week_and_Consortium,",i,"_w9w10"),
                          subset = TRUE, subset1 = paste0(i,"_w5"), subset2 = paste0(i,"_w9w10"))
      beepr::beep(3)
      print(paste0("Finished ", subs, " ", feature))
    }
  }
}

time_w5_diff_subcommunities = function(){
### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W5### - DONE
func = c("EC","KO","pathway")
  for (i in cons1[1:4]){
    for (j in cons2[1:4]){
      if (j<=i){
        next()
      }
      else{
        for (k in func[1:3]){
          feature = k
          subs = paste0(i,"_w5_vs_",j,"_w5")
          picrust_to_maaslin2( # nolint
                              variable = "Week_and_Consortium",
                              file_to_use = paste0(feature,"_merged_metagenome.tsv"),
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

time_w9w10_diff_subcommunities = function(){
  func = c("EC","KO","pathway")
  for (i in cons1[1:4]){
    for (j in cons2[1:4]){
      if (j<=i){
        next()
      }
      else{
        for (k in func[1:3]){
          feature = k
          subs = paste0(i,"_w9w10_vs_",j,"_w9w10")
          picrust_to_maaslin2( # nolint
                              variable = "Week_and_Consortium",
                              file_to_use = paste0(feature,"_merged_metagenome.tsv"),
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

same_subcommunities_diff_sex = function(){
  for (subs in cons1[1:4]){
    for (k in c("EC","KO","pathway")){
      feature = k
      picrust_to_maaslin2( # nolint
                          variable = "Sex_and_Consortium",
                          file_to_use = paste0(feature, "_merged_metagenome.tsv"),
                          feature_to_use = feature,
                          random_input = NULL,
                          reference_val = paste0("Sex_and_Consortium,Male", "_", subs),
                          subset = TRUE, subset1 = paste0("Male", "_", subs), subset2 = paste0("Female", "_", subs))
      beepr::beep(3)
      print(paste0("Finished ", subs, " ", feature))
    }
  }
}

main_function = function(){
  #subcommunities_diff_time()
  time_w5_diff_subcommunities()
  time_w9w10_diff_subcommunities()
  same_subcommunities_diff_sex()
}

main_function()
beepr::beep(8)
