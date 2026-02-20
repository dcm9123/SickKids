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
import re

# %%
def sanity_check():
    curated_cat_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_category_curation_2026.txt", sep = "\t")
    orig_pwy = curated_cat_file["Unique"].tolist()
    level1 = curated_cat_file["Level 1"].tolist()
    level1_1 = curated_cat_file["Level 1.1"].tolist()
    level1_2 = curated_cat_file["Level 1.2"].tolist()
    level1_3 = curated_cat_file["Level 1.3"].tolist()
    level2 = curated_cat_file["Level 2"].tolist()
    level2_1 = curated_cat_file["Level 2.1"].tolist()
    level2_2 = curated_cat_file["Level 2.2"].tolist()
    level2_3 = curated_cat_file["Level 2.3"].tolist()
    
    concat_list = []
    list_of_lists = [level1, level1_1, level1_2, level1_3, level2, level2_1, level2_2, level2_3]
    
    for item in list_of_lists:
        for element in item:
            element = str(element)
            concat_list.append(element)
    
    #Checked and there were no duplicates in the original file (after removing a few)
    concat_list = sorted(concat_list)
    set_concat = set(concat_list)
    set_concat = sorted(set_concat)
    
    i = 0
    j = 0
    for list in list_of_lists:
        set_list = set(list)
        if j == 0:
            i = i + 1
            j = j + 1
            print(f"Number of unique elements in the list level {i}: {len(set_list)}")
        else:
            print(f"Number of unique elements in the list level {i}_{j}: {len(set_list)}")
            j = j + 1
            if j == 4:
                j = 0
                
# %%
def annotating():
    method = ["default", "updated"]                                                 #One for each type of database used for the MinPath annotation (default MetaCyc or updated picrust2 or MetaCyc 2026)
    db = ["metacyc_2026", "picrust2_db"]                                            #One for each database used for the MinPath annotation (updated MetaCyc 2026 or PICRUSt2 database)      
    master_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/Bacterial_Metacyc_pathway_file.tsv", sep = "\t")
    master_file = master_file[["Pathways","Ontology - pathway type","Names","Common-Name-Taxa"]] #Including the bacterial master file
    master_file = master_file.drop_duplicates(subset = "Pathways", keep = "first") #Making sure that I am dropping the duplicates from the master file, so I dont get a million rows per file
    for mthd in method:     
        for database in db:
            if mthd == "default":   
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/" #Makes sure the directory is correct for the default method
                if database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                    break
            else:   #If it's not the default method, then I incorporate the two databases for the updated method
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}/"
            for file in glob.glob(os.path.join(dir1,"raw/*.txt")):      #Going over all the Minpaath reports in the raw folder of the corresponding method and database
                pwy_dictionary = {} #Initializing pwy dictionary
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
                    #print(dir1)
                os.chdir(f"{dir1}/parsed/")
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
#%%
def simplifying_categories():
    curated_cat_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_category_curation_2026.txt", sep = "\t")
    input_path = "/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/"
    os.chdir(input_path)
    method = ["updated","default"]
    db = ["metacyc_2026", "picrust2_db"]
    for mthd in method:     
        for database in db:
                if mthd == "default" and database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                    break
                elif mthd == "default" and database == db[0]:   #Makes sure the directory is correct for the default method
                    dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/parsed"
                elif mthd == "updated":
                    dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}/parsed"
                os.chdir(dir1)
                for file in glob.glob(os.path.join(dir1,"*.tsv")):
                        #k = k + 1
                        print(f"Simplifying categories for file: {file}")
                        file_in = pd.read_csv(file, sep = "\t")
                        ontology_list = file_in["Ontology - pathway type"].tolist()
                        #pwy_dict = {}
                        rows = []
                        curated_cat_set = set(curated_cat_file["Unique"].tolist())  # do this ONCE outside the loop
                        curated_cat_file_indexed = curated_cat_file.set_index("Unique")  # index by Unique for fast lookup
                        for element in ontology_list:
                            if element in curated_cat_set:
                                row = curated_cat_file_indexed.loc[element]  # get the row for this element
                                element_level1 = row["Level 1"]
                                element_level1_1 = row["Level 1.1"]
                                element_level1_2 = row["Level 1.2"]           
                                element_level1_3 = row["Level 1.3"]
                                element_level2 = row["Level 2"]
                                element_level2_1 = row["Level 2.1"]
                                element_level2_2 = row["Level 2.2"]
                                element_level2_3 = row["Level 2.3"]
                                paired_elements = list(zip([element_level1, element_level1_1, element_level1_2, element_level1_3], 
                                                           [element_level2, element_level2_1, element_level2_2, element_level2_3]))
                                sorted_pairs = sorted(paired_elements, key = lambda x: (str(x[0])=="NA", str(x[0]))) #Sorting the categories so that the NAs are always at the end of the list, and not in the middle
                                sorted_level1, sorted_level2 = zip(*sorted_pairs) #unpacks the sorted pairs into two separate lists
                                #pwy_dict[element] = {element_level1:element_level2, element_level1_1:element_level2_1, element_level1_2:element_level2_2, element_level1_3:element_level2_3}
                                #pwy_dict[element] = dict(sorted(pwy_dict[element].items(), key = lambda x: (str(x[0])=="NA", str(x[0]))))                                #level1.append(sorted([str(element_level1), str(element_level1_1), str(element_level1_2), str(element_level1_3)], key = lambda x: (x=="NA",x))) #Sorting the categories so that the NAs are always at the end of the list, and not in the middle                                    
                                #level2.append(sorted([str(element_level2), str(element_level2_1), str(element_level2_2), str(element_level2_3)], key = lambda x: (x=="NA",x)))
                            else:
                                if pd.isna(element) or element == "" or element == "NA":
                                    element_level1 = ""
                                    element_level2 = ""
                                    element_level1_1 = ""
                                    element_level1_2 = ""
                                    element_level1_3 = ""
                                    element_level2_1 = ""
                                    element_level2_2 = ""
                                    element_level2_3 = ""
                                    sorted_level1 = [element_level1, element_level1_1, element_level1_2, element_level1_3]
                                    sorted_level2 = [element_level2, element_level2_1, element_level2_2, element_level2_3]
                                    #pwy_dict[element] = {element_level1:element_level2, element_level1_1:element_level2_1, element_level1_2:element_level2_2, element_level1_3:element_level2_3}
                                    #level1.append([element_level1, element_level1_1, element_level1_2, element_level1_3])
                                    #level2.append([element_level2, element_level2_1, element_level2_2, element_level2_3])
                                else:
                                    print(f"ERROR: element {element} is not in the curated category file. Exiting now.")
                                    sys.exit()
                            rows.append({"Level 1":sorted_level1[0], "Level 1.1":sorted_level1[1], "Level 1.2":sorted_level1[2], "Level 1.3":sorted_level1[3],
                                                   "Level 2":sorted_level2[0], "Level 2.1":sorted_level2[1], "Level 2.2":sorted_level2[2], "Level 2.3":sorted_level2[3]})
                        df_out = pd.DataFrame(rows)
                        df_final = pd.concat([file_in, df_out], axis = 1)
                            #df_out = pd.DataFrame.from_dict(pwy_dict, orient = "index", columns = ["Level 1","Level 1.1","Level 1.2","Level 1.3", "Level 2","Level 2.1","Level 2.2","Level 2.3"])
                            #df_out = pd.concat(objs = [file_in, df_categories1, df_categories2], axis = 1, join = "outer", ignore_index = False)
                            #print(df_out)
                        f_out_name = file.split('/')[-1]
                            #print(f_out_name)
                            #f_out_name = re.split(r'(\d{3})', f_out_name) #To retrieve the 3 digit number in the file name, so I can keep it in the new file name and make sure I am not overwriting files
                            
                            #f_out_name = re.split(regex_val, f_out_name)[0]
                            #if database == db[0]:
                                #f_out_name = f"{f_out_name}_{mthd}_parsed_report.tsv"
                            #else:
                                #f_out_name = f"{f_out_name}_{database}_{mthd}_parsed_report.tsv"
                        print(f_out_name)    
                        df_final.to_csv(f_out_name, sep = "\t", index = False)
                        #Sanity check
                        #if k == 10:
                        #    sys.exit()

# %%
def main():
    sanity_check()
    annotating()
    simplifying_categories()

main()

# %%
