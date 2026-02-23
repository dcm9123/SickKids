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
    # Load curated category file and extract each level column as a list
    curated_cat_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_category_curation_2026.txt", sep = "\t")
    level1 = curated_cat_file["Level 1"].tolist()
    level1_1 = curated_cat_file["Level 1.1"].tolist()
    level1_2 = curated_cat_file["Level 1.2"].tolist()
    level1_3 = curated_cat_file["Level 1.3"].tolist()
    level2 = curated_cat_file["Level 2"].tolist()
    level2_1 = curated_cat_file["Level 2.1"].tolist()
    level2_2 = curated_cat_file["Level 2.2"].tolist()
    level2_3 = curated_cat_file["Level 2.3"].tolist()
    
    # Flatten all level lists into a single list for global uniqueness checking
    concat_list = []
    list_of_lists = [level1, level1_1, level1_2, level1_3, level2, level2_1, level2_2, level2_3]
    
    for item in list_of_lists:
        for element in item:
            element = str(element)
            concat_list.append(element)
    
    # Sort and deduplicate the full list (checked: no duplicates in curated file after manual curation)
    concat_list = sorted(concat_list)
    set_concat = set(concat_list)
    set_concat = sorted(set_concat)
    
    # Print unique element count per level to confirm no unexpected duplicates
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
    # Define annotation methods and databases to loop over
    method = ["default", "updated"]                                                 #One for each type of database used for the MinPath annotation (default MetaCyc or updated picrust2 or MetaCyc 2026)
    db = ["metacyc_2026", "picrust2_db"]                                            #One for each database used for the MinPath annotation (updated MetaCyc 2026 or PICRUSt2 database)      
    # Load master pathway file and keep only relevant metadata columns
    master_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/Master_metacyc_pathway_file.tsv", sep = "\t")
    master_file = master_file[["Pathways","Ontology - pathway type","Names","Common-Name-Taxa","Classification"]] #Including the bacterial master file
    #master_file = master_file.drop_duplicates(subset = "Pathways", keep = "first") #Making sure that I am dropping the duplicates from the master file, so I dont get a million rows per file
    for mthd in method:     
        for database in db:
            # Set the correct input directory based on method and database
            if mthd == "default":   
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/" #Makes sure the directory is correct for the default method
                if database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                    break
            else:   #If it's not the default method, then I incorporate the two databases for the updated method
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}/"
            # Iterate over all raw MinPath report .txt files
            for file in glob.glob(os.path.join(dir1,"raw/*.txt")):
                pwy_dictionary = {} # Initialize per-file pathway dictionary
                print(f"Processing file: {file}")
                f_out = pd.DataFrame(columns = ["Pathways","Naive","MinPath","Total families","Families found"])
                with open(file,"r") as f:
                    # Initialize counters for sanity check at the end
                    naive_yes = 0
                    naive_no = 0
                    minpath_yes = 0
                    minpath_no = 0
                    lines = f.readlines()
                    for line in lines:
                        line = line.strip()
                        pwy = line.split(" ")[1] # Extract pathway ID (2nd field)
                        # Store naive prediction: 1=found by naive, 0=not found
                        if "naive 1" in line:
                            pwy_dictionary.setdefault(pwy, {})["Naive"] = "Yes"
                        elif "naive 0" in line:
                            pwy_dictionary.setdefault(pwy, {})["Naive"] = "No"
                        # Store MinPath (parsimonious) prediction: 1=found, 0=not found
                        if "minpath 1" in line:
                            pwy_dictionary.setdefault(pwy, {})["MinPath"] = "Yes"
                        elif "minpath 0" in line:
                            pwy_dictionary.setdefault(pwy, {})["MinPath"] = "No"
                        # Store family coverage counts for each pathway
                        pwy_dictionary.setdefault(pwy, {})["Total families"] = line.split("  ")[4]
                        pwy_dictionary.setdefault(pwy, {})["Families found"] = line.split("  ")[6]
                                            
                    # Count totals across all pathways to verify data integrity
                    for value in pwy_dictionary.values():
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
                    
                    # Sanity check #1: naive total must equal minpath total (same number of pathways assessed)
                    if sum != sum2:
                        print(f"Error in file {file}: the sum of the categories does not match the total number of lines processed. Exiting now")
                        sys.exit()
                    # Sanity check #2: number of unique pathways must match total lines processed
                    if len(pwy_dictionary.keys()) != sum:
                        print(f"Error in file {file}: the number of unique pathways does not match the total number of lines processed. Exiting now")
                        sys.exit()
                # Change to the parsed output directory
                os.chdir(f"{dir1}/parsed/")
                # Convert dictionary to DataFrame and add Pathways column from index
                f_out = pd.DataFrame.from_dict(data = pwy_dictionary, orient = "index")
                f_out.reset_index(inplace = True)
                f_out.rename(columns = {"index":"Pathways"}, inplace = True)
                # Merge with master file to add pathway metadata (left join keeps all pathways)
                df_merged = f_out.merge(master_file, on = "Pathways", how = "left")
                # Build output filename based on method and database
                if database == db[0]:
                    f_out_name = f"{file.split('/')[-1].split('_minpath')[0]}_{mthd}_parsed_report.tsv"
                else:
                    f_out_name = f"{file.split('/')[-1].split('_minpath')[0]}_{database}_{mthd}_parsed_report.tsv"
                df_merged.to_csv(f_out_name, sep = "\t")
                pwy_dictionary = {} # Reset for next file
#%%
def simplifying_categories():
    # Load curated category lookup file for ontology classification
    curated_cat_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_category_curation_2026.txt", sep = "\t")
    input_path = "/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/"
    os.chdir(input_path)
    method = ["updated","default"]
    db = ["metacyc_2026", "picrust2_db"]
    # Pre-build set and indexed DataFrame ONCE for fast O(1) lookups inside the loop
    curated_cat_set = set(curated_cat_file["Unique"].tolist())
    curated_cat_file_indexed = curated_cat_file.set_index("Unique")
    for mthd in method:     
        for database in db:
                # Set input directory based on method and database combination
                if mthd == "default" and database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                    break
                elif mthd == "default" and database == db[0]:   #Makes sure the directory is correct for the default method
                    dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/parsed"
                elif mthd == "updated":
                    dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}/parsed"
                os.chdir(dir1)
                # Iterate over all parsed TSV files in the directory
                for file in glob.glob(os.path.join(dir1,"*.tsv")):
                        print(f"Simplifying categories for file: {file}")
                        file_in = pd.read_csv(file, sep = "\t")
                        ontology_list = file_in["Ontology - pathway type"].tolist()
                        rows = [] # Collect one category row dict per pathway
                        for element in ontology_list:
                            if element in curated_cat_set:
                                # Retrieve all 8 category columns in a single row lookup
                                row = curated_cat_file_indexed.loc[element]
                                element_level1 = row["Level 1"]
                                element_level1_1 = row["Level 1.1"]
                                element_level1_2 = row["Level 1.2"]           
                                element_level1_3 = row["Level 1.3"]
                                element_level2 = row["Level 2"]
                                element_level2_1 = row["Level 2.1"]
                                element_level2_2 = row["Level 2.2"]
                                element_level2_3 = row["Level 2.3"]
                                # Pair level1 and level2 so they stay aligned during sorting
                                paired_elements = list(zip([element_level1, element_level1_1, element_level1_2, element_level1_3], 
                                                           [element_level2, element_level2_1, element_level2_2, element_level2_3]))
                                # Sort alphabetically by level1, keeping "NA" entries at the end
                                sorted_pairs = sorted(paired_elements, key = lambda x: (str(x[0])=="NA", str(x[0])))
                                # Unpack sorted pairs back into two aligned level1/level2 tuples
                                sorted_level1, sorted_level2 = zip(*sorted_pairs)
                            else:
                                # Handle NaN, empty, or missing ontology entries — fill with empty strings
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
                                else:
                                    # Unexpected element not in curated file — stop and report
                                    print(f"ERROR: element {element} is not in the curated category file. Exiting now.")
                                    sys.exit()
                            # Append the sorted category values as one row dict per pathway
                            rows.append({"Level 1":sorted_level1[0], "Level 1.1":sorted_level1[1], "Level 1.2":sorted_level1[2], "Level 1.3":sorted_level1[3],
                                                   "Level 2":sorted_level2[0], "Level 2.1":sorted_level2[1], "Level 2.2":sorted_level2[2], "Level 2.3":sorted_level2[3]})
                        # Build category DataFrame and concatenate with original parsed file columns
                        df_out = pd.DataFrame(rows)
                        df_final = pd.concat([file_in, df_out], axis = 1)
                        # Write output back to the same filename in the parsed directory
                        f_out_name = file.split('/')[-1]
                        print(f_out_name)    
                        df_final.to_csv(f_out_name, sep = "\t")
                        
def counting_classification():
    # Define methods and databases to iterate over
    method = ["default", "updated"]                                                 #One for each type of database used for the MinPath annotation (default MetaCyc or updated picrust2 or MetaCyc 2026)
    db = ["metacyc_2026", "picrust2_db"]                                            #One for each database used for the MinPath annotation (updated MetaCyc 2026 or PICRUSt2 database)      
    input_path = "/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/"
    curated_file = pd.read_csv("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/metacyc_category_curation_2026.txt", sep = "\t")
    categories_1 = curated_file.loc[:,["Level 1","Level 1.1","Level 1.2","Level 1.3"]].values.tolist()
    categories_2 = curated_file.loc[:,["Level 2","Level 2.1","Level 2.2","Level 2.3"]].values.tolist()
    categories_set1 = set()
    categories_set2 = set()

    for items in categories_1:
        for element in items:
            element = str(element)
            categories_set1.add(element)
    for items in categories_2:
        for element in items:
            element = str(element)
            categories_set2.add(element)

    #for thing in sorted(categories_set1):
    #    print(thing)
    #for thing in sorted(categories_set2):
    #    print(thing)
    
    os.chdir(input_path)
    for mthd in method:     
        for database in db:
            name = []
            pwy_count = []
            naive = []
            minpath = []
            bacteria = []
            non_bacteria = []
            
            # Set input directory based on method and database
            if mthd == "default" and database == db[1]:   #Makes sure that it only runs once, as I don't have two databases for the default method
                break
            elif mthd == "default" and database == db[0]:   #Makes sure the directory is correct for the default method
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/parsed"
            elif mthd == "updated":
                dir1 = f"/Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/{mthd}/{database}/parsed"
            os.chdir(dir1)
            
            all_pwy_dict_1 = []  # list of dicts, one per file
            all_pwy_dict_2 = []  # list of dicts, one per file
            
            for file in glob.glob(os.path.join(dir1,"*.tsv")):
                pwy_dict_1 = {cat: 0 for cat in categories_set1}  # fresh count per file
                pwy_dict_2 = {cat: 0 for cat in categories_set2}  # fresh count per file
                print(f"Counting classifications for file: {file}")
                file_in = pd.read_csv(file, sep = "\t")
                name.append(file.split('/')[-1])
                pwy_count.append(len(file_in["Pathways"]))
                naive.append(len(file_in[file_in["Naive"]=="Yes"]))
                minpath.append(len(file_in[file_in["MinPath"]=="Yes"]))
                bacteria.append(len(file_in[file_in["Classification"]=="Bacteria"]))
                non_bacteria.append(len(file_in[file_in["Classification"]!="Bacteria"]))
                
                # Count level1 categories per file using column names (not fragile indices)
                for _, row in file_in.iterrows():
                    values1 = [str(row["Level 1"]), str(row["Level 1.1"]),
                                str(row["Level 1.2"]), str(row["Level 1.3"])]
                    values2 = [str(row["Level 2"]), str(row["Level 2.1"]),
                                str(row["Level 2.2"]), str(row["Level 2.3"])]
                    for category in categories_set1:
                        if category in values1:
                            pwy_dict_1[category] += 1
                    for category in categories_set2:
                        if category in values2:
                            pwy_dict_2[category] += 1
                
                # Append this file's counts to the accumulator lists
                all_pwy_dict_1.append(pwy_dict_1)
                all_pwy_dict_2.append(pwy_dict_2)
            
            # Build output AFTER all files processed — one row per file
            df_out_part1 = pd.DataFrame({"ID":name, "Total pathways":pwy_count, "Naive":naive, 
                                            "MinPath":minpath, "Bacteria":bacteria, "Non-Bacteria":non_bacteria})
            df_out_part2 = df_out_part1.copy()
            
            # Convert list of per-file dicts into DataFrames and concatenate
            df_cats1 = pd.DataFrame(all_pwy_dict_1)  # rows = files, cols = categories
            df_cats2 = pd.DataFrame(all_pwy_dict_2)
            df_out_part1 = pd.concat([df_out_part1.reset_index(drop=True), df_cats1.reset_index(drop=True)], axis=1)
            df_out_part2 = pd.concat([df_out_part2.reset_index(drop=True), df_cats2.reset_index(drop=True)], axis=1)
            
            df_out_part1.to_csv(f"Summary_classification_counting_{mthd}_{database}_level1.tsv", sep="\t", index=False)
            df_out_part2.to_csv(f"Summary_classification_counting_{mthd}_{database}_level2.tsv", sep="\t", index=False)                         
                    
# %%
def main():
    #sanity_check()             # Verify curated category file has no unexpected duplicates
    #annotating()               # Parse raw MinPath reports and merge with master pathway metadata
    #simplifying_categories()   # Add Level 1/2 category columns to each parsed file
    counting_classification()
main()

# %%
