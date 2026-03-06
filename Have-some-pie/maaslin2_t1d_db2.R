#Maaslin2
#July 21st, 2025

library(Maaslin2)
library(beepr)

picrust_to_maaslin2<-function(path_picrust2, file_to_use, variable, normalization_method,maaslin_path, 
                              transformation,metadata_file,feature_to_use,random_input,reference_val,
                              subset,subset1,subset2){
  setwd(path_picrust2)
  maaslin_input = read.table(file_to_use, header = TRUE, row.names = 1, sep='\t')
  metadata_input = read.csv(metadata_file, header = TRUE, row.names = 1)
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


#Comparing week and consortia (NS1, NS6, S2, S5) and w5, w9w10. Reference NS1_w5 # nolint

#Metadata variables: ID,Sex,Community,Subcommunity,Timepoint,Merged_weeks,Week_and_consortia,Sex_and_consortia,
#Features: EC_merged_metagenome.tsv, KO_merged_metagenome.tsv
#Normalization: TSS
#Transformation: LOG
#Random effect: Sex


### SAME SUBCOMMUNITY AND DIFFERENT TIMEPOINTS ANALYSIS### - DONE
cons1 = c("NS1","NS6","S2","S5")
cons2 = c("NS1","NS6","S2","S5")
picrust2_path_current = "/Users/danielcm/Desktop/Sickkids/PICRUSt2.3/"

for (i in cons1[1:4]) # DONE
{
  for (k in c("EC","KO"))
  {
    feature = k
    subs = paste0(i,"_w5w6_vs_",i,"_w9w10")
    picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                        variable = "Week_and_consortia",
                        normalization_method = "TSS",
                        transformation = "LOG",
                        file_to_use = paste0(feature,"_merged_metagenome.tsv"),
                        metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                        feature_to_use = feature,
                        random_input = "Sex",
                        reference_val = paste0("Week_and_consortia,",i,"_w9w10"),
                        maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                        subset = TRUE, subset1 = paste0(i,"_w5w6"), subset2 = paste0(i,"_w9w10"))
    beepr::beep(3)
    print(paste0("Finished ", subs, " ", feature))
  }
}

### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W5### - DONE
func = c("EC","KO")
for (i in cons1[1:4])
{
  for (j in cons2[1:4])
  {
    if (j<=i)
    {
      next()
    }
    else{
      for (k in func[1:2])
      {
        feature = k
        subs = paste0(i,"_w5w6_vs_",j,"_w5w6")
        picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                            variable = "Week_and_consortia",
                            normalization_method = "TSS",
                            transformation = "LOG",
                            file_to_use = paste0(feature,"_merged_metagenome.tsv"),
                            metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                            feature_to_use = feature,
                            random_input = "Sex",
                            reference_val = paste0("Week_and_consortia,",j,"_w5w6"),
                            maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/",# nolint
                            subset = TRUE, subset1 = paste0(i,"_w5w6"), subset2 = paste0(j,"_w5w6"))
        beepr::beep(3)
        print(paste0("Finished ", subs, " ", feature))
      }
    }
  }
}

### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W9W10### - DONE
cons1 = c("NS1","NS6","S2","S5")
cons2 = c("NS1","NS6","S2","S5")
func = c("EC","KO")
for (i in cons1[1:4])
{
  for (j in cons2[1:4])
  {
    if (j<=i)
    {
      next()
    }
    else{
      for (k in func[1:2])
      {
        feature = k
        subs = paste0(i,"_w9w10_vs_",j,"_w9w10")
        picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                            variable = "Week_and_consortia",
                            normalization_method = "TSS",
                            transformation = "LOG",
                            file_to_use = paste0(feature,"_merged_metagenome.tsv"),
                            metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                            feature_to_use = feature,
                            random_input = "Sex",
                            reference_val = paste0("Week_and_consortia,",j,"_w9w10"),
                            maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/",# nolint
                            subset = TRUE, subset1 = paste0(i,"_w9w10"), subset2 = paste0(j,"_w9w10"))
        beepr::beep(3)
        print(paste0("Finished ", subs, " ", feature))
      }
    }
  }
}

### SAME SUBCOMMUNITIES DIFFERENT SEX ANALYSIS ### - DONE
for (subs in c("NS1", "NS6", "S2", "S5"))
{
  for (feature in c("EC","KO"))
  {
    picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                        variable = "Sex_and_consortia",
                        normalization_method = "TSS",
                        transformation = "LOG",
                        file_to_use = paste0(feature, "_merged_metagenome.tsv"),
                        metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                        feature_to_use = feature,
                        random_input = NULL,
                        reference_val = paste0("Sex_and_consortia,Male", "_", subs),
                        maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                        subset = TRUE, subset1 = paste0("Male", "_", subs), subset2 = paste0("Female", "_", subs))
    beepr::beep(3)
    print(paste0("Finished ", subs, " ", feature))
  }
}

### NOW FOR METACYC ###
### SAME SUBCOMMUNITY AND DIFFERENT TIMEPOINTS ANALYSIS ### - DONE
    feature = "pathway"
    for (i in (cons1[1:4]))
    {
    picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                        variable = "Week_and_consortia",
                        normalization_method = "TSS",
                        transformation = "LOG",
                        file_to_use = "Pathway_merged_metagenome.tsv",
                        metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                        feature_to_use = "pathway", # Change to pathway when doing metacyc
                        random_input = "Sex",
                        reference_val = paste0("Week_and_consortia,",i,"_w9w10"),
                        maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                        subset = TRUE, subset1 = paste0(i,"_w5w6"), subset2 = paste0(i,"_w9w10"))
    beepr::beep(3)
    print(paste0("Finished ", i, " ", feature))
    }

### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W5W6 ### - DONE
for (i in cons1[1:4])
{
  for (j in cons2[1:4])
  {
    if (j<=i)
    {
      next()
    }
    else{
      feature = k
      subs = paste0(i,"_w5w6_vs_",j,"_w5")
      picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                          variable = "Week_and_consortia",
                          normalization_method = "TSS",
                          transformation = "LOG",
                          file_to_use = "Pathway_merged_metagenome.tsv",
                          metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                          feature_to_use = "pathway", # Change to pathway when doing metacyc
                          random_input = "Sex",
                          reference_val = paste0("Week_and_consortia,",j,"_w5w6"),
                          maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                          subset = TRUE, subset1 = paste0(i,"_w5w6"), subset2 = paste0(j,"_w5w6"))
      beepr::beep(3)
      print(paste0("Finished ", i," vs ",j, " ", feature))
      }
  }
}

### DIFFERENT SUBCOMMUNITIES AT THE SAME TIMEPOINT W9W10 ### - DONE
for (i in cons1[1:4])
{
  for (j in cons2[1:4])
  {
    if (j<=i)
    {
      next()
    }
    else{
      feature = "pathway"
      subs = paste0(i,"_w9w10_vs_",j,"_w9w10")
      picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                          variable = "Week_and_consortia",
                          normalization_method = "TSS",
                          transformation = "LOG",
                          file_to_use = "Pathway_merged_metagenome.tsv",
                          metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                          feature_to_use = "pathway", # Change to pathway when doing metacyc
                          random_input = "Sex",
                          reference_val = paste0("Week_and_consortia,",j,"_w9w10"),
                          maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                          subset = TRUE, subset1 = paste0(i,"_w9w10"), subset2 = paste0(j,"_w9w10"))
      beepr::beep(3)
      print(paste0("Finished ", i," vs ",j," ", feature))
      }
  }
}

### SAME SUBCOMMUNITY AND DIFFERENT SEX ANALYSYS ###
for (i in cons1[1:4]){
  feature = "pathway"
  subs = i
    picrust_to_maaslin2(path_picrust2 = picrust2_path_current, # nolint
                        variable = "Sex_and_consortia",
                        normalization_method = "TSS",
                        transformation = "LOG",
                        file_to_use = "Pathway_merged_metagenome.tsv",
                        metadata_file = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/maaslin2_july2025/metadata_ps_without_w7.csv", # nolint
                        feature_to_use = feature,
                        random_input = NULL,
                        reference_val = paste0("Sex_and_consortia,Male", "_", subs),
                        maaslin_path = "/Users/danielcm/Desktop/Sickkids/Maaslin2.4/", # nolint
                        subset = TRUE, subset1 = paste0("Male", "_", subs), subset2 = paste0("Female", "_", subs))
    beepr::beep(3)
    print(paste0("Finished ", subs, " ", feature))
  }
