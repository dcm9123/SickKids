# Daniel Castaneda Mogollon, PhD
# February 18th, 2026
# Purpose: This script was made to verify the categories of the pathways in the MetaCyc 2026 database.
# by simply comparing the number of lines that have a key word from the parent categories in Metacyc

# Expected:
# Bioluminescence - 13
# Biosynthesis - 2186
# Degradation - 1157
# Detoxification - 70
# Generation of Precursor Metabolites and Energy - 265
# Gylcan pathways - 188
# Macromolecule modification - 82
# Metabolic clusters - 72
# Metabolic activation/inactivation/interconversion - 53
# Signal transduction pathways - 0
# Superpathways - 383
# Transport - 1


# SCRIPT DID NOT WORK

#%%
category_file = open("/Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/Pathway_categories_metacyc_2026.txt", "r")
biol = []
biosyn = []
degr = []
detox = []
gen_prec = []
glycan = []
macromod = []
metabolic_clusters = []
metabolic_act = []
signal = []
superpwy = []
transport = []

with category_file as f:
    for line in f:
        line = line.strip()
        line = line.lower()
        if "super-pathways" in line: #good
            superpwy.append(line)
        elif "bioluminescence" in line: # good 
            biol.append(line)
        elif "glycan" in line: #good
            glycan.append(line)
        elif "detox" in line:
            detox.append(line)
        elif "secondary-biosynthesis" in line or "energy-" in line:
            gen_prec.append(line)
        elif "biosynthesis" in line:
            biosyn.append(line)
        elif "degradation" in line:
            degr.append(line)
        elif "macromolecule" in line:
            macromod.append(line)
        elif "metabolic clusters" in line:
            metabolic_clusters.append(line)
        elif "activation-inactivation" in line:
            metabolic_act.append(line)
        elif "signal transduction" in line: #good
            signal.append(line)
        elif "transport" in line: #good
            transport.append(line)

print(f"Bioluminescence - {len(biol)}")
print(f"Biosynthesis - {len(biosyn)}")
print(f"Degradation - {len(degr)}")
print(f"Detoxification - {len(detox)}")
print(f"Generation of Precursor Metabolites and Energy - {len(gen_prec)}")
print(f"Glycan pathways - {len(glycan)}")
print(f"Macromolecule modification - {len(macromod)}")
print(f"Metabolic clusters - {len(metabolic_clusters)}")
print(f"Metabolic activation/inactivation/interconversion - {len(metabolic_act)}")
print(f"Signal transduction pathways - {len(signal)}")
print(f"Superpathways - {len(superpwy)}")
print(f"Transport - {len(transport)}")
# %%
