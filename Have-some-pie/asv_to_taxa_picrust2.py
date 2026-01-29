# Daniel Castaneda Mogollon, PhD
# September 25th, 2024
# Purpose: To link the PICRUST2 output table containing feature names and ASV contribution per sample
# to the taxonomy table from FemMicro16S and ran through Phyloseq to get the taxonomy of each ASV.

import pandas as pd
import os
import argparse
from time import sleep

parser = argparse.ArgumentParser(description='Link PICRUST2 output features/functions with taxonomy table.', add_help=True)
parser.add_argument('--feature_type', type=str,choices=['KO','EC','pathway'],help='Type of feature to analyze: KO, EC, or pathway', required=True) #nargs passes a list back to the code with + being one or more values
parser.add_argument('--features', type=str, nargs='+', required=True,help='List of KOs, ECs, or pathways of interest separated by spaces, e.g., K00001 K00002; EC1.2.3.4 EC2.3.1.4') #nargs passes a list back to the code with + being one or more values')
parser.add_argument('--main_attribute', type=str, choices=['sex','consortia','week'], required=True,
                    help='Main attribute to work with, e.g., week for (w5w6 vs w9w10), or sex (male vs female),or consortia.'
                    ' Must be one of the three options exactly as shown here.')
parser.add_argument('--pairwise_comparison', required=True,choices = ['sex','consortia','week'],help='Provide the attribute that you '
                    'will be comparing the main attribute against, i.e. week (if you chose consortia as main attribute, '
                    'then the comparison will be done by sex in each consortia)')
# Example usage: python asv_to_taxa_picrust2.py --KOs K00001 K00231 K00345 --main_attribute consortia --pairwise_comparison week

def revieweing_arguments(): #This method reviews the arguments passed by the user, and ensures they are valid.
    args = parser.parse_args()
    list_features = args.features 
    main_attribute = args.main_attribute
    pairwise_comparison = args.pairwise_comparison
    function_type = args.feature_type #This will be either KO, EC, or pathway. I only need the first value since the user should only provide one type of function at a time.

    if main_attribute == "consortia": #Defines the list that will be used based on the main attribute, i.e. if 'consortia', it means that each consortia will be the main component where a pairwise comparison will be made.
        list_to_work_with = ['ns1','ns6','s2','s5']
    elif main_attribute == "week":
        list_to_work_with = ['w5','w9w10']
    elif main_attribute == "sex":
        list_to_work_with = ['male','female']
    else:
        print("Error: Invalid main_attribute argument. Please choose 'sex', 'consortia', or 'week'. Terminating program.")
        exit()

    if(main_attribute == pairwise_comparison): #The same attribute cannot be used for both main and pairwise comparison, i.e. cannot compare consortia as main and as pairwise.
        print("Error: main_attribute and pairwise_comparison cannot be the same. Please choose different attributes. Terminating program.")
        exit()
    
    print("\n")
    print("Main attribute:", main_attribute)
    print("Function type to analyze:", function_type)
    print("Features to analyze:", list_features)
    print("Pairwise comparison attribute:", pairwise_comparison)
    print("\n")

    sleep(10) #Sleep for 5 seconds to give the user time to read the arguments provided
    
    return(list_to_work_with, list_features, function_type, pairwise_comparison)

def taxa_name_harmonization_and_grouping(list_to_work_with,list_features,function_type):
    path = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions" #Path to the folder containing the output folders from PICRUSt2
    df_dict = {}
    for element in list_to_work_with: #This means I should have folders for all possible comparisons, i.e. ns1_output, w5w6_output, male_output, etc.
        print(f"Processing element: {element}")
        stratified_output = pd.read_csv(f"{path}/{element}_output/{element}_{function_type}_metagenome_out/pred_metagenome_strat.tsv", sep="\t") #Reading the stratified output from PICRUst2
        taxa_table = pd.read_csv(f"{path}/{element}_output/ps_{element}_asv_final_renamed_tax_table.csv", sep=",") #Reading the taxonomy table from Phyloseq output, remember to delete the "" in bash and change the first column name to 'sequence'
        df_merged = pd.merge(left = stratified_output,right = taxa_table,on="sequence",how="left")  #Merge both of the data frames based on the sequence (which is the ASV ID)
        df_merged['species_final'] = df_merged['species_final'].fillna('sp.') #Make a new column that fills in the missing species with 'sp.'
        df_merged['genus_and_species'] = df_merged['genus_final'] + " " + df_merged['species_final'] #Make a new column that joins the genus and species
        #df_merged['taxa_merging'] = np.where(df_merged['species_final'] == 'sp.', df_merged['genus_final'], df_merged['genus_and_species']) #Make a new column that merges the genus and species but if species is missing, just use genus
        for feature in list_features: #Check which KOs are in the list provided by the user
            if feature not in df_merged['function'].values:
                print(f"Feature {feature} not found in the file. Skipping.")
                continue    #Iterates to the next value of the for loop
            else:
                print(f"Processing feature: {feature}")
                feature_df = df_merged[df_merged['function'] == feature] #Subset the merged dataframe to only include rows with the KO of interest
                feature_df = feature_df.groupby(['genus_and_species', 'function'], as_index=False).sum()
                feature_comb = feature+"_" + element
                df_dict[feature_comb] = {'dataframe': feature_df, 'element': element, 'feature':feature} #Store the dataframe in a dictionary with the feature as key
        print("\n")
                
    return(df_dict)

def normalizing_contribution_by_taxa(df_dict, pairwise_comparison):
    for key,item in df_dict.items():
        data_frame = item['dataframe']
        element = item['element']
        feature = item['feature']
        print(key,item)
        if pairwise_comparison == 'week':
            path = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2/Picrust2_predictions" #Path to the folder containing the output folders from PICRUSt2        if pairwise_comparison == "week":
            dfw5w6 = data_frame.iloc[:,0:2] #Subset the dataframe to first include the first three columns (genus_and_species, function and sequence)
            dfw9w10 = data_frame.iloc[:,0:2] #Subset the dataframe to first include the first three columns (genus_and_species, function and sequence)
            w5w6_columns = data_frame.filter(regex='week5').columns
            w9columns = data_frame.filter(regex='week9').columns
            w10_columns = data_frame.filter(regex='week10').columns
            dfw5w6 = pd.concat([dfw5w6,data_frame[w5w6_columns]], axis=1)
            dfw9w10 = pd.concat([dfw9w10,data_frame[w9columns],data_frame[w10_columns]], axis=1)
            dfw5w6['Abundance'] = dfw5w6.iloc[:,2:].sum(axis=1) #Sum the abundance across all samples for that particular KO
            dfw9w10['Abundance'] = dfw9w10.iloc[:,2:].sum(axis=1) #Sum the abundance across all samples for that particular KO
            dfw5w6['Proportion'] = (dfw5w6['Abundance'] / dfw5w6['Abundance'].sum()) * 100
            dfw9w10['Proportion'] = (dfw9w10['Abundance'] / dfw9w10['Abundance'].sum()) * 100
            dfw5w6.to_csv(f"{path}/{element}_output/{feature}_taxa_table_contribution_by_{pairwise_comparison}_w5w6.csv", index=False)
            dfw9w10.to_csv(f"{path}/{element}_output/{feature}_taxa_table_contribution_by_{pairwise_comparison}_w9w10.csv", index=False)
            
            #dfw5w6.to_csv(f"{path}/)
                
                # feature_df[f"Bacterial_contribution_in_{element}"] = (feature_df['abundance'] / feature_df['abundance'].sum()) * 100
                
                #feature_df.to_csv(f"{path}/{element}_output/{feature}_taxa_table_contribution_by_{pairwise_comparison}.csv", index=False)
                #print(f"Output saved to {path}/{element}_output/KO_of_interest_taxa_table_contribution.csv")


#The next step is to divide each .csv file by their KO, and add a new argument that asks the user what kind of comparison to make.
#That comparison will allow me to create a new column that has status A vs B (i.e. week5 vs week9). I then can normalize the percentage
#of contribution of each taxa for that particular KO.

arguments = revieweing_arguments()
new_df_list = taxa_name_harmonization_and_grouping(arguments[0], arguments[1],arguments[2])
normalizing_contribution_by_taxa(new_df_list, arguments[3])


