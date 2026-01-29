# Daniel Castaneda Mogollon, PhD
# September 30th, 2025
# Purpose: This script was created to summarize the results of the Maaslin2 analysis.

import glob
import pandas as pd
import os

df_out = pd.DataFrame()
path = '/Users/danielcm/Desktop/diammatics/T1D/Maaslin2/'
os.chdir(path) 

df_list = []
for dir in glob.glob("*/"):
    dir_path = os.path.join(path, dir)
    os.chdir(dir_path)                  #Takes me to each individual directory 
    results=[]
    for dir2 in glob.glob("*/"):
        dir2_path = os.path.join(dir_path, dir2)
        os.chdir(dir2_path)             #Takes me to each individual directory within each category of Maaslin2 by group comparison
        for file in glob.glob("*all_results.tsv"):
            if dir=="sex/":
                random_effect = None
            else:
                random_effect = "Sex"
            variable = dir.replace('/','')  # Extract variable name from directory name 
            value1 = dir2.split('_')[5:7]  # Extract comparison name from directory name
            value2 = dir2.split('_')[8:10]
            comparison = ''.join(value1) + ' vs ' + ''.join(value2)
            feature = dir2.split('_')[4]  # Extract feature type from directory name
            df = pd.read_csv(file, sep='\t')
            #random_effect = dir2.split('')
            significant_results1 = len(df[(df['qval'] < 0.05) & ((df['coef'] < -1) | (df['coef'] > 1))]) #Gets the number of significant results based on qval and coef thresholds
            significant_results2 = len(df[(df['qval'] < 0.001) & ((df['coef'] < -1) | (df['coef'] > 1))]) #Gets the number of significant results based on qval and coef thresholds 
            sign_left = len(df[(df['qval'] < 0.05) & (df['coef'] < -1)])
            sign_right = len(df[(df['qval'] < 0.05) & (df['coef'] > 1)])
            sign_but_not_enrich= len(df[(df['qval']<0.05) & ((df['coef'] > -1) & (df['coef'] < 1))]) #Gets the number of significant results based on qval threshold but not enriched by log2FC
            not_significant = len(df[df['qval'] >= 0.05]) #Gets the number of non-significant results based on qval threshold
            mice_number = int(df['N'].mean()) #Gets the number of mice in the comparison, ignore the mean function.
            total_features = df['feature'].nunique()
            results.append({'Feature':feature, 
                            'Random effect':random_effect, 
                            'Fixed effect':variable, 
                            'Comparison':comparison, 
                            'Significant results (qval<0.05 & |coef|>1)':significant_results1, 
                            'Significant results (qval<0.001 & |coef|>1)':significant_results2, 
                            'Enriched left ' +' (qval<0.05 & coef<-1)':sign_left, 
                            'Enriched right ' +' (qval<0.05 & coef>1)':sign_right,
                            'Significant but not enriched (qval<0.05 & |coef|<1)':sign_but_not_enrich,
                            'Not significant (qval>=0.05)':not_significant, 
                            'Total features':total_features, 
                            'Mice number':mice_number})
    df_out = pd.DataFrame(results)
    df_out.to_csv(f'/Users/danielcm/Desktop/diammatics/T1D/Maaslin2/{dir}{variable}_enrichment_feature_summary.tsv', sep='\t', index=False)
print(len(results))
            
            
