# Daniel Castaneda Mogollon, PhD
# January 29th, 2026
# Script designed to merge multiple MetaCyc pathway, reaction, and enzyme files into a single file for easier processing

import pandas as pd
import os
import re

def filtering_by_bacteria():
    
    os.chdir("/Users/danielcm/Desktop/diammatics/T1D/MetaCyc_files/")
    df = pd.read_csv("pathways_with_species.txt", sep="\t")
    df['Common-Name'] = df['Common-Name'].str.strip()
    df['Common-Name2'] = df['Common-Name2'].str.strip()

    df['Common-Name'] = df['Common-Name'].str.replace(' // ',',')
    df['Common-Name2'] = df['Common-Name2'].str.replace(' // ',',')
    
    
    regex_criteria = re.compile(r'(?i)(\bBacteria\b|<bacteria>)|\b(?:Pseudomonadota|Actinomycetota|Bacillota|Cyanobacteriota|Gammaproteobacteria|Enterobacteriaceae|Bacteroidales|Lactobacillales|Streptomycetaceae|Micromonosporales|Kitasatosporales)\b|\b(?:Escherichia|Salmonella|Pseudomonas|Streptomyces|Bacillus|Clostridium|Enterococcus|Klebsiella|Acinetobacter|Rhodococcus|Cupriavidus|Mycobacterium|Corynebacterium|Vibrio|Yersinia|Azotobacter|Bifidobacterium|Lacticaseibacillus|Lactiplantibacillus|Phocaeicola)\s+[a-z][a-z-]+(?:\s+\w+)?\b', re.IGNORECASE)
    regex_criteria2 = re.compile(r'\b(?:Archaea|Eukaryota|Viridiplantae|Embryophyta|Magnoliopsida|Metazoa|Mammalia|Chordata|Vertebrata|Fungi|Ascomycota|Basidiomycota|Opisthokonta|Chlorophyta|Bryophyta|Stramenopiles|Alveolata|Amoebozoa|Apicomplexa|Euglenozoa|Haptophyta|Bacillariophyta|Oomycota)\b', re.IGNORECASE)
    
    # Add na=False to all .str.contains() calls
    bacteria_present = df['Common-Name'].str.contains(regex_criteria, na=False) | df['Common-Name2'].str.contains(regex_criteria, na=False)
    no_bacteria = df['Common-Name'].str.contains(regex_criteria2, na=False) | df['Common-Name2'].str.contains(regex_criteria2, na=False)
     
    df_bacteria = df[bacteria_present]
    df_non_bacteria = df[~bacteria_present & no_bacteria]
    df_neither = df[~bacteria_present & ~no_bacteria]
    
    print(df_bacteria.shape)
    print(df_non_bacteria.shape)
    print(df_neither.shape)
    #print(df_neither['Common-Name2'])
    print(df.shape)
    return df

#def reading_files():

def global_function():
    df = filtering_by_bacteria()
    print(df['Common-Name'].head())
    
global_function()
