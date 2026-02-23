# Daniel Castaneda Mogollon, PhD
# February 2nd, 2026
# This script will take the *paid* database from MetaCyc and prepare the files with a more complete framework for Minpath
# MinPath requires a file with two columns where the 1st one is the pathway ID and the 2nd one is the enzyme ID reaction (EC.2.3.4.1)

#%%
import pandas as pd
import os


#%%
path = "/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/"
os.chdir(path)

#%%
pathway_file = pd.read_csv("Master_Metacyc_pathway_file.tsv", sep = "\t")

#%%
def cleaning_ECs(df):
    df = df.copy() # to avoid modifying the original dataframe
    print(f"The number of retrieved Pathways from MetaCyc's 2026 file is: {df.shape[0]}")
    df = df[df['Classification'] == "Bacteria"]
    print(f"The original number of Pathways for bacteria is: {df.shape[0]}")
    df = df.dropna(subset=["EC-Number"])
    print(f"The number of Pathways with an existing EC number(s) is: {df.shape[0]}")
    for col in ["EC-Number","Reaction-List"]:
        if col == "EC-Number":
            df[col] = df[col].str.replace("EC-","")
        df[col] = df[col].astype(str)
        df[col] = df[col].str.split(" // ")
        #df[col] = df[col].tolist()
    #display(df)
    return(df)

def getting_rxns(df):
    df = df.copy() # to avoid modifying the original dataframe
    print(f"The number of retrieved Pathways from MetaCyc's 2026 file is: {df.shape[0]}")
    df = df[df['Classification'] == 'Bacteria']
    print(f"The original number of Pathways for bacteria is: {df.shape[0]}")
    df = df.dropna(subset=["Reactions of pathway"])
    print(f"The number of Pathways with an existing reaction list is: {df.shape[0]}")
    df['Reactions of pathway'] = df.apply(lambda row: row['Reactions of pathway'].split(" // "), axis=1)
    return(df)

#%%
def exploding_pairs(df):
    tmp = df.copy()
    #tmp['Paired-values'] = tmp.apply(lambda row: list(zip(row['EC-Number'], row['Reaction-List'])), axis=1)
    tmp = tmp.explode("Reactions of pathway", ignore_index = True)
    #tmp[["EC-Number","Reaction-List"]] = pd.DataFrame(tmp["Paired-values"].tolist(), index=tmp.index)
    #tmp = tmp.drop(columns=["Paired-values"])
    tmp.to_csv("Bacterial_Metacyc_pathway_file.tsv", sep = "\t", index = False)
    return(tmp)

def refining_for_minpath(df):
    df_minpath = df[["Pathways","Reactions of pathway"]].copy()
    df_minpath.to_csv("Minpath_ready_Metacyc_pathway_file.tsv", sep = "\t", index = False, header = False)
    return(df_minpath)

#%%
def main():
    #df = cleaning_ECs(pathway_file)
    df2 = getting_rxns(pathway_file)
    df3 = exploding_pairs(df2)
    df_final = refining_for_minpath(df3)
    display(df_final)
    
if __name__ == "__main__":
    main()

# %%
