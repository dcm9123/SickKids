# Daniel Castaneda Mogollon, PhD
# February 2nd, 2026
# Script to prepare MinPath input files from the annotated EggNOGs (ECs) to get the metabolic pathways This script takes the EC count from all of my genomes
# That file was generated before for PICRUSt2. Then, it makes an output file where the first colummn is the genome name with a random ID identifier (i.e. S_NS1_Bf_084_1, S_NS1_Bf_004_2 ...
# and the right column has the E.C. identified for that genome. In many many instances there are multiple EC counts per genome, so the output repeats the name of the genome with a new identifier
# each time, and the exact same EC ID. (i.e) This format is necessary for MinPath to run with my local file. No headers are needed for this to work and no 'EC:' identifier is required.
# S_NS1_Bf_084_1    1.14.19.3
# S_NS1_Bf_084_2    1.14.19.3
# S_NS1_Bf_084_3    1.23.4

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
