# Daniel Castaneda Mogollon, PhD
# 14:19
# eggnog_sanity_check.py.py
# Purpose: This script compares the results of the merged file vs the individual EggNOG files to make sure
# all the numbers are matching properly.

import os
import pandas as pd
import glob
import re

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)  # Let pandas auto-detect width

path = "/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/picrust2_june232025/eggnogs_annotations/"    #path where the eggnog emapper.annotation files are
os.chdir(path)
dir1 = "run1/"                  #directory of first run
dir2 = "run2/"                  #directory of 2nd run
dir3 = "NCBI_flagged/"          #directory of NCBI genomes that were fixed
dir_final = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/"            #Change this when needed! directory where the final merged KO and EC files are located
df_output = "EggNOG_sanity_check.tsv" #Output file for the entire workflow

KO_merged_final_file = pd.read_csv(dir_final+"KO_for_picrust2_renamed.tsv",sep="\t")  #Read the file I used for picrust2 input (the final one)
KO_merged_final_file["assembly"] = KO_merged_final_file["assembly"].apply(lambda x: x.split("_")[0:4]) #Separates the four instances of the sample ID name (removes v1v9,etc)
KO_merged_final_file["assembly"] = KO_merged_final_file["assembly"].apply(lambda x: "_".join(x)) #Joines the four elements that were separated

EC_merged_final_file = pd.read_csv(dir_final+"EC_for_picrust2_renamed.tsv",sep="\t")  #Same but for ECs
EC_merged_final_file["assembly"] = EC_merged_final_file["assembly"].apply(lambda x: x.split("_")[0:4])
EC_merged_final_file["assembly"] = EC_merged_final_file["assembly"].apply(lambda x: "_".join(x))

file_dict = {}                                      #dictionary for file name and values of EC, and KOs
ec_identifier = re.compile(r'\d+\.\d+\.\d+\.\d+')   #This are regex for counting ECs with 4 digits and 3 dots
ko_identifier = re.compile(r'ko:K\d{5}')            #This regex counts ko:K and five digits after
dir_list = [dir1,dir2,dir3]                         #Iterating over the three directories where the emapper.annotations are
for dir in dir_list:
    for file in glob.glob(dir+"*.annotations"):
        file_m = re.split("/", file)[1] #Takes away the run1/ run2/ NCBI_flagged/ part of the string
        file_m = re.split("\.", file_m)[0]  #Renaming for the dictionary key name by taking away the .emapper.annotations
        with open(file,"r") as f:
            count_ko = 0
            count_ec = 0
            lines = f.readlines()
            for line in lines:
                count_ko = len(re.findall(ko_identifier, line)) + count_ko
                count_ec = len(re.findall(ec_identifier, line)) + count_ec
            file_dict[file_m] = count_ko,count_ec                                  #Writes the counts of KO and ECs per sample from the individual files

df_ko_count_individual = pd.DataFrame.from_dict(file_dict,orient='index',columns=['KO','EC'])

ko_sums = KO_merged_final_file.set_index("assembly").apply(pd.to_numeric, errors='coerce').sum(axis=1)
ko_sums = pd.DataFrame(ko_sums, columns=["KO merged count"])
ec_sums = EC_merged_final_file.set_index("assembly").apply(pd.to_numeric, errors='coerce').sum(axis=1)
ec_sums = pd.DataFrame(ec_sums, columns=["EC merged count"])

ko_sums.reset_index(inplace=True)
ec_sums.reset_index(inplace=True)

df_final = pd.DataFrame.from_dict(file_dict, orient='index', columns=["KO original count", "EC original count"])
df_final.index.name = "Assembly"
df_final.reset_index(inplace=True)

df_confirmation = df_final.merge(ko_sums, left_on="Assembly", right_on="assembly").merge(ec_sums, left_on="Assembly", right_on="assembly")
df_confirmation.drop(columns=["assembly_x","assembly_y"],inplace=True)
matches_ko = df_confirmation["KO original count"] == df_confirmation["KO merged count"]
matches_ec = df_confirmation["EC original count"] == df_confirmation["EC merged count"]
df_confirmation["KO match"],df_confirmation["EC match"] = matches_ko,matches_ec

print("Files for sanity check written to: /Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions") #Change this if needed!
ko_sums.to_csv("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/ko_merged_count.tsv",sep="\t",index=False) #This too!
ec_sums.to_csv("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/ec_merged_count.tsv",sep='\t',index=False) #And this!
df_final.to_csv("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/individual_sums.tsv",sep="\t",index=False)   #And this!
df_confirmation.to_csv("/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions/confirmation_matching.tsv",sep="\t",index=False) #Final output file with the comparison
print("There are a few ECs (not KOs) that do not match. This is because I need to fix the regex in the individual file counting for ECs. Check the eggnog_sanity_check.xlsx file for details.")