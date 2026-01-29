### Daniel Castaneda Mogollon, PhD
### October 1st, 2025
### Purpose: This script was created to fetch the MetaCyc pathway names and classification from the individual output files from Maaslin2.
### output. 

#Run this in the conda environment 'selenium' where I have installed the selenium package. Type python directly, don't
#use the default interpreter from VSCode

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import json, time
import os
import pandas as pd
import requests
from bs4 import BeautifulSoup
import glob

def fetch_pathway_info(pathway_id):
    url = f"http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=PATHWAY&object={pathway_id}"
    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36"} # Mimic a real browser request
    response = requests.get(url, headers=headers) # Fetch the page content
    soup = BeautifulSoup(response.text, "html.parser") # Parse the HTML content
    name = soup.title.text.strip() if soup.title else "" # Get the title as the pathway name
    classification = "" # Initialize classification
    for p in soup.find_all("p", class_="ecoparagraph"): # Find the paragraph with classification info
        if "Superclasses:" in p.text: # Look for the right paragraph
            tds = p.find_all("td") # Get all table data cells
            if len(tds) > 1: # Ensure there is a second cell
                classes = [a.text.strip() for a in tds[1].find_all("a")] # Extract class names
                classification = " > ".join(classes) # Join class names with ' > '
            break
    return {"PathwayID": pathway_id, "Name": name, "Classification": classification} # Return the collected info as a dictionary

file_to_fill = "/Users/danielcm/Desktop/diammatics/T1D/PICRUSt2.2/Pathway_merged_metagenome.tsv"
df = pd.read_csv(file_to_fill, sep="\t", index_col=0, header=0)
df = pd.DataFrame(df)
df = df.transpose()
categories = {}
f_out = open("metacyc_pathway_details2.tsv", "w")
f_out.write("PathwayID\tName\tClassification\n")
#for i,pathway in enumerate(df.index, start = 1): # Uncomment this line if you want to do it on a file where the pathways are the columns
pathway_list = ["THRESYN-PWY","TRNA-CHARGING-PWY","TRPSYN-PWY","UBISYN-PWY","UDPNAGSYN-PWY","VALSYN-PWY",
                "P281-PWY","PWY-3781","CRNFORCAT-PWY","P162-PWY","PROPFERM-PWY","PWY-2201","PWY-6148"]
for i, pathway in enumerate(pathway_list, start=1):
    print(f"Processing pathway: {pathway}")
    data = fetch_pathway_info(pathway)
    categories[pathway] = data
    f_out.write(f"{data['PathwayID']}\t{data['Name']}\t{data['Classification']}\n")
    delay = 3 + (i % 4)
    print(f"Waiting for {delay} seconds to avoid rate limiting...")
    time.sleep(delay)
    if i % 20 == 0:
        print("Taking a longer break of 8 seconds...")
        time.sleep(8)
    print(categories)
f_out.close()
print("Done")

'''
base_path = "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2" # Base directory containing subdirectories
for path1 in glob.glob(os.path.join(base_path, "*/")): # Iterate over each subdirectory in the base path
    if path1 == "/Users/danielcm/Desktop/diammatics/T1D/Maaslin2/sex/":
        print("Skipping sex folder")
        continue
    for path2 in glob.glob(os.path.join(path1, "*_pathway_*/")):
        work_path = os.path.abspath(path2)
        print("Processing file in:", work_path)
        all_results_file = os.path.join(work_path, "all_results.tsv")
        if not os.path.exists(all_results_file):
            print(f"File not found: {all_results_file}")
            continue
        df = pd.read_csv(all_results_file, sep="\t")
        pathway_ids = df['feature'].str.replace(".", "-") # Replace dots with hyphens in pathway IDs
        results = []
        for i, pathway in enumerate(pathway_ids, 1):
            data = fetch_pathway_info(pathway)
            results.append(data)
            print(f"Processing {i}: {pathway} out of {len(pathway_ids)}. Pathway classification: {data.get('Classification', '')}")
            delay = 3 + (i % 4)
            print(f"Waiting for {delay} seconds to avoid rate limiting...")
            time.sleep(delay)
            if i % 20 == 0:
                print("Taking a longer break of 8 seconds...")
                time.sleep(8)
        results_df = pd.DataFrame(results)
        classification_cols = results_df['Classification'].str.split(' > ', expand=True)
        classification_cols.columns = [f'Class{i+1}' for i in range(classification_cols.shape[1])]
        results_df = pd.concat([results_df, classification_cols], axis=1)
        results_df.to_csv(os.path.join(work_path, "metacyc_pathway_info.tsv"), sep="\t", index=False)
        '''