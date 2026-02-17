# Daniel Castaneda Mogollon, PhD
# February 5th, 2026
# Purpose: This script was made to parse each MinPath report and detailed report into one single
# file containing the info from all of the genomes analyzed.
# Input: The directories where the MinPath reports are located
# Output: A .tsv file with all the genomes and their predicted pathways by naïve approach, parsimonious, and the 
# default and updated MetaCyc databases.

#%%
import pandas as pd
import os
import glob as glob
import sys

# %%
method = ["default", "updated"]                                                 #One for each type of database used for the MinPath annotation (default MetaCyc or updated picrust2 or MetaCyc 2026)
db = ["metacyc_2026", "picrust2_db"]                                            #One for each database used for the MinPath annotation (updated MetaCyc 2026 or PICRUSt2 database)      
master_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/Bacterial_Metacyc_pathway_file.tsv", sep = "\t")
master_file = master_file[["Pathways","Ontology - pathway type","Names","Common-Name-Taxa"]] #Including the bacterial master file
master_file = master_file.drop_duplicates(subset = "Pathways", keep = "first") #Making sure that I am dropping the duplicates from the master file, so I dont get a million rows per file
# %%
for mthd in method:     
    for database in db:
        if mthd == "default":   
            dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}" #Makes sure the directory is correct for the default method
            if database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                break
        else:   #If it's not the default method, then I incorporate the two databases for the updated method
            dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}"
        pwy_dictionary = {} #Initializing pwy dictionary
        i = 0
        for file in glob.glob(os.path.join(dir1,"raw/*.txt")):      #Going over all the Minpaath reports in the raw folder of the corresponding method and database
            #print(f"Processing file: {file}\n")
            print(f"Processing file: {file}")
            f_out = pd.DataFrame(columns = ["Pathways","Naive","MinPath","Total families","Families found"])    #Naming columns of interest
            with open(file,"r") as f:
                naive_yes = 0
                naive_no = 0
                minpath_yes = 0
                minpath_no = 0
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    pwy = line.split(" ")[1] #Split the text file so I retrieve the pwy id
                    if "naive 1" in line:
                        pwy_dictionary.setdefault(pwy, {})["Naive"] = "Yes"
                    elif "naive 0" in line:
                        pwy_dictionary.setdefault(pwy, {})["Naive"] = "No"
                    if "minpath 1" in line:
                        pwy_dictionary.setdefault(pwy, {})["MinPath"] = "Yes"
                    elif "minpath 0" in line:
                        pwy_dictionary.setdefault(pwy, {})["MinPath"] = "No"
                    pwy_dictionary.setdefault(pwy, {})["Total families"] = line.split("  ")[4]
                    pwy_dictionary.setdefault(pwy, {})["Families found"] = line.split("  ")[6]
                                        
                for value in pwy_dictionary.values(): #Getting which pwys are present by each method of Minpath (naive vs parsimonious) 
                    if value["Naive"] == "Yes":
                        naive_yes = naive_yes + 1
                    elif value["Naive"] == "No":
                        naive_no = naive_no + 1
                    if value["MinPath"] == "Yes":
                        minpath_yes = minpath_yes + 1
                    elif value["MinPath"] == "No":
                        minpath_no = minpath_no + 1
                    sum = naive_yes + naive_no
                    sum2 = minpath_yes + minpath_no
                
                #Sanity check #1
                if sum != sum2:         #Counting them to make sure that the sum of the categories matches the total number of lines processed and the number of unique pathways
                    print(f"Error in file {file}: the sum of the categories does not match the total number of lines processed. Exiting now")
                    sys.exit()
                #Sanity check #2
                if len(pwy_dictionary.keys()) != sum:
                    print(f"Error in file {file}: the number of unique pathways does not match the total number of lines processed. Exiting now")
                    sys.exit()
            #Changing directory to the parsed dir
            os.chdir(f"{dir1}/parsed")
            f_out = pd.DataFrame.from_dict(data = pwy_dictionary, orient = "index")
            f_out.reset_index(inplace = True)
            f_out.rename(columns = {"index":"Pathways"}, inplace = True)
            df_merged = f_out.merge(master_file, on = "Pathways", how = "left") #Merging master file with each individual parsed file
            if database == db[0]:
                f_out_name = f"{file.split('/')[-1].split('_minpath')[0]}_{mthd}_parsed_report.tsv"
            else:
                f_out_name = f"{file.split('/')[-1].split('_minpath')[0]}_{database}_{mthd}_parsed_report.tsv"
            df_merged.to_csv(f_out_name, sep = "\t")
            pwy_dictionary = {} 
        #break  
# %%
