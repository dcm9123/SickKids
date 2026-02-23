# Daniel Castaneda Mogollon, PhD
# February 12th, 2026
# Purpose: This script will take the EC numbers from the annotated genomes and conver them to reactions
# by using PICRUSt2 file of 'ec_level4_to_metacyc_reactions' (or something like that).
# Input: A pathway where all of the annotated genomes are located with their EC numbers.
# Input2: The PICRUSt2 file with the EC to reaction mapping.
# Output: A file for each genome with the reactions per corresponding EC number.


#%%
import pandas as pd
import os
import glob as glob
import sys

#%%
path = "/Users/danielcm/Desktop/SickKids/MetaCyc/"
os.chdir(path)

picrust2_db = pd.read_csv("Master_Files/picrust2_ec_level4_to_metacyc_rxn.tsv", sep = "\t", header=None)
picrust2_db[0] = picrust2_db[0].apply(lambda row: row.replace("EC:",""))
methods = ["Sanger","Barrnap"]
list_ecs_picrust2 = picrust2_db[0].tolist()


for method in methods:
    for file in glob.glob("EC_annotated_genomes/ECs/"+method+"/*.tsv"):                 #Going over all EC-annotated files by both methods
        ecs_not_present = []
        print(f"Processing file: {file}")
        new_name = file.split("/")[-1].replace("ecs","reactions")                       #Making a new name for the output .tsv and report files
        f_out = open(path+"EC_annotated_genomes/Reactions/"+method+"/"+new_name, "w")   #Opening output file .tsv to write the reactions
        f_out_report = open(path+"EC_annotated_genomes/Reactions/"
                            +method+"/"+new_name.replace(".tsv","_report.txt"), "w")    #Opening output file .txt to write the ECs that are not present in the PICRUSt2 database
        print(f"Processing file {file}")
        df_input = pd.read_csv(file, sep = "\t", header=None)
        for row in df_input.itertuples(index = False, name = None):                     #Going over all rows of the EC-annotated input file
            ec_input = row[1]                                                           #Getting the EC number from the second column of the input file
            if (picrust2_db.iloc[:,0] == ec_input).any():                               #If picrust2 db has that EC across all the rows of its first column, then...
                reaction_list = picrust2_db.loc[picrust2_db[0]==ec_input,1].to_list()   #Save the reactions into a list of one element (so far)...
                if("," in reaction_list[0]):                                            #If there are commas, then we can assume there are multiple reactions in this one-element list
                    reactions = reaction_list[0].split(",")                             #Divide the 'one-element' list into a proper list with multiple reactions
                    for reaction in reactions:
                        f_out.write(str(row[0])+"\t"+reaction+"\n")                     #If there are multiple reactions, then write one row for each with the same ID of the original row of the input file
                else:
                    f_out.write(str(row[0])+"\t"+reaction_list[0]+"\n")                 #If not, then just write one row with the original ID and the only reaction in the list
            else:
                ecs_not_present.append(ec_input)
                f_out_report.write(str(ec_input) + "\n")
        print(f"There are a total of {len(ecs_not_present)} ECs not present in the PICRUSt2 database from file {file}")
        f_out.close()
        f_out_report.close()
            

        # %%
