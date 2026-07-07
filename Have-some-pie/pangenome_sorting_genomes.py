# Daniel Castaneda Mogollon, PhD
# July 7th, 2026
# Purpose: This script takes the input files from the pangenome fasta files (core, accessory, and strain-specific)
# and returns a sorted output file by genome ID

from email import header
import re
import os
import glob as glob

path = "/Users/danielcm/Desktop/SickKids/Anvio/07_PANGENOMICS/"
os.chdir(path)

for dir in glob.glob("*PROJECT"):
    os.chdir(path+dir)
    print(f"Processing directory: {dir}")
    for file in glob.glob("*_genes.fa"):
        with open(file, "r") as f:
            dictionary = {}
            filename = os.path.basename(file)
            directory = os.path.dirname(file)
            print(f"Processing file: {filename} in directory: {directory}")
            lines = f.readlines()
            for line in lines:
                if line.startswith(">"):
                    if "genome_name:" not in line:
                        print(f"Warning: skipping malformed header in {filename}: {line.strip()}")
                        exit(1)
                    header = line.strip()
                    line = header.split("genome_name:")[1]
                    genome_id = line.split("|")[0]
                    dictionary[header] = [genome_id, ""]
                else:
                    dictionary[header][1] = dictionary[header][1] + line.strip()
        dictionary = dict(sorted(dictionary.items(), key = lambda x: x[1][0]))
        
        f_out_name = filename.replace(".fa","_sorted.fa")
        with open(f_out_name, "w") as f_out:
            for key, value in dictionary.items():
                f_out.write(f"{key}\n")
                f_out.write(f"{value[1]}\n")

print("Sorting complete.")


#print(dictionary.keys())
