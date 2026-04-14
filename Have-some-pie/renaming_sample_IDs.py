# Daniel Castaneda Mogollon, PhD
# 04/13/2026
# Purpose: To rename the old IDs given from the Danska Lab with the SRA IDs
# submitted to NCBI. Then use the new count table with the new IDs to run PICRUSt2.

#%%
import os
import pandas as pd

#%%

path = "/Users/danielcm/Desktop/SickKids/"
os.chdir(path)

comm = "S5" #Change this for the community you want to work with
metadata_file_granato = pd.read_csv("Metadata/Danska_diabetes_metadata364_20260409.csv", sep=',', header=0)
file_to_read = f"Phyloseq2/Filtered_to_use/{comm}_filtered_ASVs_count600_len400_prev20_f.csv" #This is the formatted ASV table
file_to_read2 = f"Phyloseq2/{comm}_filtered_ASVs_count600_len400_prev20.csv" #In here I use the seq and ID
file_to_write = f"Phyloseq2/Filtered_to_use/{comm}_filtered_ASVs_count600_len400_prev20_f_SRA.csv" #Output file

#%%
my_dict = {} #This dictionary is meant to replace the 'plate1_ns1_w5_R_L_S23' with the SRA ID
asv_dict = {} #This dictionary is meant to replace the 'p2_asv13' with the actual sequence
df_in = pd.read_csv(file_to_read, sep=',', header=0)
df_in2 = pd.read_csv(file_to_read2, sep=',', header=0)

for row in metadata_file_granato.itertuples():
    my_dict[row.STATA_FemMicro364] = row.SRA_sample_name

for row in df_in2.itertuples(): #Unlock for NS1, NS6
    asv_dict[row.asv_id] = row.asv_seq

df_in.rename(columns=my_dict, inplace=True)
df_in["asv_id"] = df_in["asv_id"].map(asv_dict) #Unlock for NS1, NS6
df_in.to_csv(file_to_write, sep=',', index=False)


#%%