# Daniel Castaneda Mogollon, PhD
# February 2nd, 2026
# This script will take the *paid* database from MetaCyc and prepare the files with a more complete framework for Minpath
# MinPath requires a file with two columns where the 1st one is the pathway ID and the 2nd one is the enzyme ID reaction (EC.2.3.4.1)

#%%
import pandas as pd
import os


#%%
path = "/Users/danielcm/Desktop/SickKids/"
os.chdir(path)

#%%
pathway_file = pd.read_csv("../diammatics/T1D/Metacyc_files/Master_Metacyc_pathway_file.tsv", sep = "\t")

#%%
def cleaning_ECs(df):
    df = df.copy() # to avoid modifying the original dataframe
    print(f"The number of retrieved Pathways from MetaCyc's 2025 file is: {df.shape[0]}")
    df = df[df['Classification'] == "Bacteria"]
    print(f"The original number of Pathways for bacteria is: {df.shape[0]}")
    df = df.dropna(subset=["EC-Number"])
    print(f"The number of Pathways with an existing EC number is: {df.shape[0]}")
    for col in ["EC-Number","Reaction-List"]:
        if col == "EC-Number":
            df[col] = df[col].str.replace("EC-","")
        df[col] = df[col].astype(str)
        df[col] = df[col].str.split(";")
        #df[col] = df[col].tolist()
    #display(df)
    return(df)

#%%
def exploding_pairs(df):
    tmp = df.copy()
    tmp['Paired-values'] = tmp.apply(lambda row: list(zip(row['EC-Number'], row['Reaction-List'])), axis=1)
    tmp = tmp.explode("Paired-values", ignore_index = True)
    tmp[["EC-Number","Reaction-List"]] = pd.DataFrame(tmp["Paired-values"].tolist(), index=tmp.index)
    tmp = tmp.drop(columns=["Paired-values"])
    tmp.to_csv("Exploded_Metacyc_pathway_file.tsv", sep = "\t", index = False)
    return(tmp)

def refining_for_minpath(df):
    df_minpath = df[["Pathways","EC-Number"]]
    df_minpath.to_csv("Minpath_ready_Metacyc_pathway_file.tsv", sep = "\t", index = False, header = False)
    return(df_minpath)

#%%
def main():
    df = cleaning_ECs(pathway_file)
    df2 = exploding_pairs(df)
    df_final = refining_for_minpath(df2)
    display(df_final)
    
if __name__ == "__main__":
    main()

# %%
