# Daniel Castaneda Mogollon, PhD
# Purpose: To merge the MetaCyc pathway information with the Maaslin2 results for better interpretation.
# Date: October 2th, 2025

import pandas as pd
import glob
import os
import numpy as np

path = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2"
for path1 in glob.glob(os.path.join(path, "*/")):
    for working_path in glob.glob(os.path.join(path1, "*_pathway_*/")):
        os.chdir(working_path)
        df_maaslin = pd.read_csv("all_results.tsv", sep="\t")
        df_metacyc = pd.read_csv("metacyc_pathway_info.tsv", sep="\t")
        df_maaslin['feature'] = df_maaslin['feature'].str.replace(".", "-", regex=False)
        df_maaslin.rename(columns={'feature': 'PathwayID'}, inplace=True)
        df_merged = pd.merge(df_maaslin, df_metacyc, on="PathwayID",how='left')
        if 'NegLog10(qval)' not in df_merged.columns:
            df_merged['NegLog10(qval)'] = -1 * np.log10(df_merged['qval'])
            
        if df_merged.shape[0] != df_maaslin.shape[0]:
            print(f"Warning: Mismatch in number of rows after merging in {working_path}")
        else:
            df_merged.to_csv("all_results_with_metacyc_info.tsv", sep="\t", index=False)
            print(f"Successfully merged data!")

exit()

    