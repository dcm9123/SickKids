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

i = 1

for method in methods:
    for file in glob.glob("EC_annotated_genomes/ECs/"+method+"/*.tsv"):
        i = i+1
        print(f"Processing file: {file}")
        new_name = file.split("/")[-1].replace("ecs","reactions")
        #print(f"Processing file {file}")
        df_input = pd.read_csv(file, sep = "\t", header=None)
        with open(path+"EC_annotated_genomes/Reactions/"+method+"/"+new_name, "w") as f_out:
            for row in df_input.itertuples(index=False, name=None):
                genome_id, ec_input = row[0], row[1]
                if (picrust2_db.iloc[:,0] == ec_input).any():      #if picrust2 db has that EC across all the rows of its first column, then...
                    reaction_list = picrust2_db.loc[picrust2_db[0]==ec_input,1].to_list()
                    if not reaction_list:
                        continue
                    for reaction_entry in reaction_list:
                        if "," in reaction_entry:
                            reactions = reaction_entry.split(",")
                            for reaction in reactions:
                                f_out.write(str(genome_id)+"\t"+reaction+"\n")
                        else:
                            f_out.write(str(genome_id)+"\t"+reaction_entry+"\n")
        if i > 2:
            sys.exit()

            

        # %%
