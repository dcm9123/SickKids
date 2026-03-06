#Daniel CM, PhD
#March 1st, 2026

#%%
import pandas as pd
import os
import glob as glob

#%%
path = "/Users/danielcm/Desktop/Sickkids/"
os.chdir(path)
# %%
subdirs = ["consortium","sex","w5","w9w10"]
master_f = "MetaCyc/Master_Files/Master_Metacyc_pathway_file_with_categories.tsv"
df_master = pd.read_csv(master_f, sep = "\t")
df_master = df_master.loc[:,["Pathways","Names","Level 1","Level 1.1","Level 2.1","Level 1.2","Level 2.2","Level 1.3","Level 2.3"]]

for file in glob.glob(f"Maaslin2.4/*/*_pathway_*/all_results.tsv"):
    df = pd.read_csv(file, sep = "\t")
    df['feature'] = df['feature'].apply(lambda x: x.replace(".","-"))
    df_merged = df.merge(df_master,how = "left", left_on="feature",right_on="Pathways")
    f_out = file.replace(".tsv","with_categories.tsv")
    df_merged.to_csv(f_out,sep = "\t")
    
    

# %%
