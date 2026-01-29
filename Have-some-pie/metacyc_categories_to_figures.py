# Daniel Castaneda Mogollon, PhD
# January 16th, 2026
# Script to map MetaCyc pathway categories to the MetaCyc IDs from heatmaps or networks

import pandas as pd
import glob

def laura_function():
    file1 = "/Users/danielcm/Desktop/diammatics/T1D/metacyc_pathway_details.tsv"
    file2 = "/Users/danielcm/Downloads/Pathways-from-All-Pathways-of-MetaCyc.txt"

    # Read the MetaCyc pathway details file
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")

    df1_IDs = df1['PathwayID'].tolist()
    df2_IDs = df2['Pathways'].tolist()
    present = []
    absent = []

    for idx, row in df2.iterrows():
        pathway_id = row['Pathways']
        if pathway_id in df1_IDs:
            present.append(1)
        else:
            present.append(0)

    for idx, row in df1.iterrows():
        pathway_id = row['PathwayID']
        if pathway_id not in df2_IDs:
            absent.append(pathway_id)

    df2['Flagged_T1D'] = present
    x = df2['Flagged_T1D'].sum() # 311 pathways overlapped
    y = df1.shape[0] # 317 pathways in T1D file

    print(f"There are a total of {x} pathways from the T1D merged metagenome in the MetaCyc all pathways file")
    print(f"{y} pathways are not present, those are: {absent}")

    for file in glob.glob("/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Heatmaps/*when_pairwise*.tsv"):
        sign_pathways = []
        found = []
        not_found = {}
        print("Working on file:", file)
        df3 = pd.read_csv(file, sep="\t", index_col=0)
        for pathway in df3.iloc[:,0]:
            sign_pathways.append(pathway)
        for ids, row in df2.iterrows():
            pathway_id = row['Pathways']
            if pathway_id in sign_pathways:
                found.append(1)
            else:
                not_found[file] = pathway_id
                found.append(0)
        df2['Flagged'+file.split("/")[-1]] = found
        print(not_found.values())                

    df2_to_csv = df2.to_csv(sep="\t", index=False)
    f_out = open("Pathways-from-All-Pathways-of-MetaCyc_dcm.tsv", "w") 
    f_out.write(df2_to_csv)
    f_out.close()
    return()

def network_categories():
    file1 = "/Users/danielcm/Desktop/diammatics/T1D/metacyc_pathway_details2.tsv"
    consortia = ["NS1", "S2"]
    timepoints = ["w5", "w9w10"]
    df1 = pd.read_csv(file1, sep="\t")
    pathway_dict = {}
    for idx, row in df1.iterrows():
        pathway_id = row['PathwayID']
        name = row['Name']
        level1 = row['Level1']
        level2 = row['Level2']
        level3 = row['Level3']
        level4 = row['Level4']
        level5 = row['Level5']
        level6 = row['Level6']
        pathway_dict[pathway_id] = {
            "Name": name,
            "level1": level1,
            "level2": level2,
            "level3": level3,
            "level4": level4,
            "level5": level5,
            "level6": level6
        }

    for consortium in consortia:
        for file in glob.glob(f"/Users/danielcm/Desktop/diammatics/T1D/Maaslin2.3/Networks/{consortium}/*network_pathway_IDs.txt"):
            print("Working on file:", file)
            df_in = pd.read_csv(file, header=None, sep = "\t")
            df_in['Pathway_ID'] = df_in[0].astype(str)
            column_names = ["Name","Level1","Level2","Level3","Level4","Level5","Level6"]
            df_in[column_names] = df_in["Pathway_ID"].map(pathway_dict).apply(pd.Series)
            f_out = file.replace("_network_pathway_IDs.txt", "_network_pathway_IDs_with_categories.tsv")
            df_in.to_csv(f_out, sep="\t", index=False)
            
            
            #NS1_w9w10_when_pairwise_is_NS1_w5_and_NS1_w9w10_network_pathway_IDs.txt

network_categories()
    