# Daniel Castaneda Mogollon, PhD
# February 2nd, 2026
# Script to prepare MinPath input files from the annotated EggNOGs (ECs) to get the metabolic pathways

#%%
import pandas as pd
import os
import glob as glob

#%%
path = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/"
os.chdir(path)

df = pd.read_csv("EC_for_picrust2_renamed.tsv", sep = "\t")

# %%
df.rename(columns = lambda col: col.replace("EC:",""), inplace=True)
# %%
i = 1
for row in df.itertuples():
    print(f"Processing genome {row.assembly}")
    f_out = open(f"EC_annotated_genomes/{row.assembly}_minpath_ecs.tsv", "w")
    i = i +1
    j = 1
    for column in range(1,len(df.columns),1):
        number = df.at[row.Index, df.columns[column]]
        number = int(number)
        print(number)
        if number > 0:
            for k in range(1,number+1,1):
                f_out.write(f"{row.assembly}_{j}\t{df.columns[column]}\n")
                j = j+1
        else:
            continue
    f_out.close()
        
    

# %%
