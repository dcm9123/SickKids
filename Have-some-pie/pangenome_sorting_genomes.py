# Daniel Castaneda Mogollon, PhD
# July 7th, 2026
# Purpose: This script takes the input files from the pangenome fasta files (core, accessory, and strain-specific)
# and returns a sorted output file by genome ID

import os
import glob as glob
import sqlite3

path = "/Users/danielcm/Desktop/SickKids/Anvio/07_PANGENOMICS/EBOLTEAE-PROJECT/"
os.chdir(path)

conn = sqlite3.connect("EBOLTEAE-PAN.db")
cur = conn.cursor()

def sort_fasta_by_genome_id():
    os.chdir(path)
    for dir in glob.glob("*PROJECT"):
        os.chdir(path+dir)
        print(f"Processing directory: {dir}")
        for file in glob.glob("*.fa"):
            if file.endswith("_sorted.fa"):
                continue  # Skip already sorted files
            else:
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
                            line = header.split("genome_name:")[1] #This parses Et_015|gene_callers_id:4355
                            genome_id = line.split("|")[0] #This parses Et_015
                            gc_id = header.split("|")[1].split(":")[1] #This parses GC_00000021
                            dictionary[header] = [genome_id, "", gc_id]
                        else:
                            dictionary[header][1] = dictionary[header][1] + line.strip()
                dictionary = dict(sorted(dictionary.items(), key = lambda x: (x[1][0], x[1][2])))
                
                f_out_name = filename.replace(".fa","_sorted.fa")
                with open(f_out_name, "w") as f_out:
                    for key, value in dictionary.items():
                        f_out.write(f"{key}\n")
                        f_out.write(f"{value[1]}\n")

    print("Sorting complete.")
    return()

def metadata_layer_tag():
    os.chdir(path)
    for dir in glob.glob("*PROJECT"):
        os.chdir(path+dir)
        dictionary = {}
        #print(f"Accessory genome file not found, proceeding with core and strain files")
        print(f"Processing directory: {dir}")
        core_file = "core_sorted.fa"
        strain_file = "strain_specific_sorted.fa"
        if os.path.exists("accessory_sorted.fa"):
            print(f"Accessory genome file found")
            accessory_file = "accessory_sorted.fa"
            files = [core_file, accessory_file,strain_file]
        else:
            print(f"Accessory genome file not found, proceeding with core and strain files")
            files = [core_file,strain_file]
        f_out_misc = "pangenome_level.txt"
        for file in files:
            with open(file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith(">"):
                        id = line.split("|")[1].split(":")[1]
                        if file=="core_sorted.fa":
                            dictionary[id] = "Core"
                        elif file=="strain_specific_sorted.fa":
                            dictionary[id] = "Strain-specific"
                        else:
                            dictionary[id] = "Accessory"
        print(dictionary)
        with open(f_out_misc, "w") as f_out:
            f_out.write("gene_cluster_id\tpangenome_level\n")
            dictionary = dict(sorted(dictionary.items(), key = lambda x: (x[1], x[0])))
            print(dictionary)
            for key, value in dictionary.items():
                f_out.write(f"{key}\t{value}\n")

    print("Metadata layer tagging complete.")

def get_genome_set_for_gc(gc_id):
    cur.execute("SELECT DISTINCT genome_name FROM gene_clusters WHERE gene_cluster_id=?", (gc_id,))
    return tuple(sorted(r[0] for r in cur.fetchall()))
#sort_fasta_by_genome_id()
#metadata_layer_tag()

def get_genome_for_gc(gc_id):
    cur.execute("SELECT DISTINCT genome_name FROM gene_clusters WHERE gene_cluster_id=?", (gc_id,))
    rows = cur.fetchall()
    return rows[0][0] if len(rows) == 1 else None  # should be exactly 1 for strain-specific

category_order = {"Core": 0, "Accessory": 1, "Strain-specific": 2}

with open("pangenome_level.txt") as f:
    next(f)
    rows = [line.strip().split("\t") for line in f]

def sort_key(row):
    gc_id, category = row
    if category == "Strain-specific":
        genome = get_genome_for_gc(gc_id)
        return (category_order[category], 0, (genome,), gc_id)
    elif category == "Accessory":
        genome_set = get_genome_set_for_gc(gc_id)
        return (category_order[category], -len(genome_set), genome_set, gc_id)
    return (category_order[category], 0, (), gc_id)

rows.sort(key=sort_key)

with open("items_order.txt", "w") as f_out:
    for gc_id, category in rows:
        f_out.write(gc_id + "\n")


