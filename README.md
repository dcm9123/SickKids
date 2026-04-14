# T1D Project Entries from 2026

### February 2nd, 2026

I have decided to create a new repository that only contains the publicly available files as well as the necessary ones (as long as they are below the capacity of Github's storage). I was asked to do an annotation profile of each genome with the ECs and therefore the MetaCyc pathways. After dealing so long with MetaCyc, Kevin and I figured out that the best way to do it is to export the 'smart tables' by indicating which fields we are interested in. However, the website is crap and we can only do 5 fields at a time, then I came up with a script `pathway_metacyc_merger.py` that takes the fields I need across multiple files, and then they are merged into one large one, called `Master_Metacyc_pathway_file.tsv`. I had to come up with a list of bacterial names, genus, family, etc in order to keep a database exclusively for bacteria, as their fields are really bad and not filtered by Kingdom.

I have finished coding the `metacyc_genome_annotator.py` which retrieves the genome number and prints the EC multiple times (n times of CNV found for that cell). The input file is the `EC_for_picrust2_renamted.tsv` and it generates an output file for every row, printing the number and ID of each EC. I made sure all rows had a different ID by simply adding a counter. Now that I have both files, MinPath should run smoothly. But that's a task for tomorrow!

### February 3rd, 2026

I ran the `MinPath` software across one sample to test it out. This is the command I used: `python MinPath.py -ec ../../diammatics/T1D/PICRUSt2.2/EC_annotated_genomes/S_NS1_Af_002_minpath_ecs.tsv -report test_report -details test_details` inside this directory: `/Users/danielcm/Desktop/SickKids/MinPath_2026`

That gave me the following output files:
a) test_details:
```
path PWY-5143 fam0 1 fam-found 1 # long-chain fatty acid activation
   6.2.1.3 hits 3 # long-chain-fatty-acid-CoA ligase
path PWY-5083 fam0 3 fam-found 1 # NAD(P)/NADPH interconversion
   2.7.1.23 hits 1 # NAD kinase
path PWY-7723 fam0 7 fam-found 1 # bacterial bioluminescence
   1.5.1.38 hits 1 # FMN reductase (NADPH)
path PWY-40 fam0 2 fam-found 2 # putrescine biosynthesis I
   4.1.1.19 hits 1 # arginine decarboxylase
   3.5.3.11 hits 1 # agmatinase
path PWY-6562 fam0 4 fam-found 2 # norspermidine biosynthesis
   1.2.1.11 hits 1 # aspartate-semialdehyde dehydrogenase
   4.1.1.96 hits 1 # carboxynorspermidine decarboxylase
path PWY-8072 fam0 1 fam-found 1 # alanine racemization
   5.1.1.1 hits 1 # alanine racemase
path COMPLETE-ARO-PWY fam0 18 fam-found 9 # superpathway of aromatic amino acid biosynthesis
   4.2.3.4 hits 1 # 3-dehydroquinate synthase
   4.2.1.10 hits 1 # 3-dehydroquinate dehydratase
   1.1.1.25 hits 1 # shikimate dehydrogenase
   2.7.1.71 hits 1 # shikimate kinase
   2.5.1.19 hits 1 # 3-phosphoshikimate 1-carboxyvinyltransferase
   4.2.3.5 hits 1 # chorismate synthase
   5.4.99.5 hits 1 # chorismate mutase
   4.2.1.51 hits 1 # prephenate dehydratase
   1.3.1.12 hits 1 # prephenate dehydrogenase
```

b) test_report:
```
path PWY-5143 any n/a  naive 1  minpath 1  fam0  1  fam-found  1  name  long-chain fatty acid activation
path PWY-5083 any n/a  naive 1  minpath 1  fam0  3  fam-found  1  name  NAD(P)/NADPH interconversion
path PWY-5114 any n/a  naive 1  minpath 0  fam0  8  fam-found  2  name  UDP-sugars interconversion
path PWY-7723 any n/a  naive 1  minpath 1  fam0  7  fam-found  1  name  bacterial bioluminescence
path PWY-40 any n/a  naive 1  minpath 1  fam0  2  fam-found  2  name  putrescine biosynthesis I
path PWY-43 any n/a  naive 1  minpath 0  fam0  3  fam-found  1  name  putrescine biosynthesis II
path PWY-6305 any n/a  naive 1  minpath 0  fam0  3  fam-found  1  name  superpathway of putrescine biosynthesis
path P101-PWY any n/a  naive 1  minpath 0  fam0  4  fam-found  1  name  ectoine biosynthesis
```


According to the developer, these files can be read like this:
```
4.3 How to read the MinPath report file?
    e.g, demo.ko.minpath
    ...
    path 00030 kegg n/a  naive 1  minpath 1  fam0  42  fam-found  18  name  Pentose phosphate pathway
    path 00031 kegg n/a  naive 1  minpath 0  fam0  12  fam-found  2  name  Inositol metabolism
    ...
    1) path 00030, path 00031 are the KEGG pathway IDs (if your input file has fig families, the pathways will then be SEED subsystems)
    2) kegg n/a, indicates pathway reconstruction of your input dataset if not available (note: this information is only available for the genomes annotated 
	 in KEGG database); for input with fig families, KEGG is replaced by SEED 
    3) naive 1 or 0: the pathway is reconstructed, or not, by the naive mapping approach
    4) minpath 1 or 0: the pathway is kept, or removed by MinPath
    5) fam0: the total number of families involved in the corresponding pathway
    6) fam-found: the total number of involved families that are annotated
    7) name: the description of the corresponding pathway (subsystem)
  
 4.4 How to read the MinPath detailed report file?
    This report file lists all the pathways found MinPath, and a list of families each pathway includes
    e.g. demo.ko.minpath.details
    ...
    path 00010 fam0 56 fam-found 27 # Glycolysis / Gluconeogenesis
       K00001 hits 6 # E1.1.1.1, adh
       K00002 hits 1 # E1.1.1.4, adh
    ...
    1) path 00010 is the KEGG pathway ID, and fam0 and fam-found are the same as in 4.3
    2) this pathway includes families, K00001 and K00002, and so on 
       (here K numbers are used for KEGG families, and FIG ids for FIG families)
    3) family K00001 has 6 hits (i.e., 6 proteins/reads are annotated as this family)
```

In summary, the `test_details` file contains only the MinPath MetaCyc pathways using the minimum parsimonious approach. In my case, if I look at COMPLETE-ARO-PWY, there are 18 families involved in that pathway (families as Kegg families, Enzymes, etc (in my case each family is the EC ID), and it found 9 of them present in there, and there is one hit of each enzyme/read found and annotated as this family. 

The `test_report` contains other type of information. Here, I also have the MetaCyc Pathway ID, a naive binary number of 0 or 1, where 0 means it is absent using a naive approach, and 1 saying it is present using the same approach. The minpath 1 or 0 is the same but using the minimum parsimonious approach. The fam0 and fam numbers are the same, and the name is the description of the corresponding pathway. I can use this file to compare the annotation differences between genomes. 

That command was *specifically* to run my ECs using the developers database! This is the modified version of it: `python MinPath.py -any ../../diammatics/T1D/PICRUSt2.2/EC_annotated_genomes/S_NS1_Af_002_minpath_ecs.tsv -map ../Minpath_ready_Metacyc_pathway_file.tsv -report test_report2 -details test_details2` inside this directory: `/Users/danielcm/Desktop/SickKids/MinPath_2026`

The difference is the mapping file. If I use the default one, I have a total of 5,009 and mine has 7,187 pathways (with the bacteria filter and the latest db available from MetaCyc). 
When using both commands, I have the following results:

default pathways found from MetaCyc: 211
My db pathways found: 420

Here is the difference between the naïve and the parsimonious approach in a figure *that does not belong to me*:
https://www.google.com/url?sa=t&source=web&rct=j&url=https%3A%2F%2Fjournals.plos.org%2Fploscompbiol%2Farticle%3Fid%3D10.1371%2Fjournal.pcbi.1000465&ved=0CBYQjRxqFwoTCPDrnfSyvpIDFQAAAAAdAAAAABBO&opi=89978449<img width="600" height="407" alt="image" src="https://github.com/user-attachments/assets/433f6a0f-9d34-4e39-9573-a0ccbfc8ab38" />

That image shows in a nutshell the differences between these two approaches. MinPath will use all of the families of ECs present and find the path that explains the lowest number of metabolic pathways possible, whereas the naïve approach overestimates this number by saying a pathway is present if a family is present for that pathway.

### February 4th, 2026
After running the following command `for file in EC_annotated_genomes/*; do python MinPath.py -any ${file} -map Minpath_ready_Metacyc_pathway_file.tsv -report ${file%%.tsv}_report.txt -details ${file%%.tsv}_detailed_report.txt; done;` to get all of the results and details across all of the genomes, I put the files in the directory `/Users/danielcm/Desktop/SickKids/MinPath_2026/annotated_pathways_genomes/`. I also made new directories to simplify my results by type of report. The next step is to add a summary and the category of the pathway next to each finding. I also need to parse the data from the report into a .tsv format to make it more 'friendly' to the user. 

### February 5th, 2026
Today I modified my scripts so I could have a better database for whomever wants to use it. I have uploaded the scripts into MetaCyc of my SickKids repository. So far, I have done the following:


1. Download the MetaCyc smart tables by pathways, reactions, enzymes (this was hard... and slow... but I got it done).
2. Then I wrote a script called `metacyc_genome_annotator.py` that takes the ECs from each genome (annotated by Prokka and EggNOG) and retrieves the genome number and the number of ECs and ID associated with that genome. I did this across all genomes.
3. Then I wrote a script called `pathway_metacyc_merger.py` that merges the multiple files downloaded from MetaCyc into one large file containing only the fields I need (Pathway ID, Enzymes, Reactions, Categories, Summary, amongs many others). This file is called `Master_Metacyc_pathway_file.tsv` and it has all of the metabolic pathways from MetaCyc, not just bacteria. In the same script, I added a function that retrieves the bacterial pathways (and perhaps some eukaryotic ones that happen to also encode for that bacterial pathway). This file is called `Master_Metacyc_pathway_file.tsv`. 
4. Then I wrote a script called `prepping_minpath_files.py` that uses the `Master_Metacyc_pathway_file.tsv` to generate a data frame with only bacterial pathways and the associated enzymes (ECs). This file has two columns, one with the pathway ID and the other with ECs associated with that pathway. The difference is that this file does not need headers, and it is expects to have multiple identical pathway IDs if there are multiple ECs associated with that pathway (i.e. if a pathway has 3 ECs, then this pathway ID is found 3 times in the left column and each EC next to it.) This is the file that MinPath needs to run properly and it is called `Bacterial_Metacyc_pathway_file.tsv`. This is needed for MinPath to run properly.
5. Finally, I ran MinPath across all genomes using the command mentioned yesterday.

Today I am going to work on categorizing the pathways into broader categories (e.g., Metabolism, Genetic Information Processing, Environmental Information Processing, etc.) based on the MetaCyc classification. This will help in understanding the functional capabilities of each genome at a higher level. I will also try to parse the MinPath output files into a more user-friendly format, possibly converting them into .tsv files for easier analysis and visualization.

### February 10th, 2026 

Working on the Journal Club paper "Gut microbial production of imidazole propionate drives Parkinson's pathologies" (2025) by Huyunji Park et al, nature communications.

### February 11th, 2026

I'll start working on the script that aims to subset the PICRUSt2's MetaCyc database. I will be annotating my genomes based on these approaches:
1) Default MetaCyc database in MinPath
2) The 'Master' MetaCyc database I created from the 2026 download.
3) The PICRUSt2's MetaCyc database that I will subset based on the bacterial pathways I have in my 'Master' MetaCyc database.

### February 12th, 2026

I realized some of the pathways I downloaded and present in the Master MetaCyc file do not have the corresponding ECs in them, which is a problem because MinPath needs the ECs to run. I will run a code to check which of these are missing. My suspicion is that most of them belong to Superpathways, as the corresponding 'ECs' are listed as subpathways, so I will have to remove the 'ECs' and replace them with reaction IDs which are present in all the pathways of my Master file table.

### February 13th, 2026

I have finalized a code called `picrust2_metacyc_annotation.py` that takes the input annotated EC files for each Barrnap and Sanger-annotated genomes, and it retrieves the ECs so it can write the matching reaction taken from the PICRUSt2 embedded database. If there are multiple reactions, it writes one row for each reaction with the same ID of the original row of the input file. If there are ECs that are not present in the PICRUSt2 database, it writes them in a separate report file and prints the total number of ECs not present in the PICRUSt2 database for each input file. This way, I can keep track of which ECs are missing and how many reactions I have in total. Naturally, all of my reaction files should be larger in rows than my EC files. It will look like this:

EC-input file:

|ID|EC|
|:--- | :--- |
| S_NS1_Af_002_v1v9_1 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_2 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_3 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_4 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_5 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_6 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_7 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_8 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_9 | 	1.6.5.3 |
| S_NS1_Af_002_v1v9_10 |	1.6.5.3 |
| S_NS1_Af_002_v1v9_11 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_12 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_13 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_14 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_15 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_16 |	3.6.4.12 |
| S_NS1_Af_002_v1v9_17 |	3.2.1.23 |

New output reaction file:

|ID|Reaction|
|:--- | :--- |
| S_NS1_Af_002_v1v9_1 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_2 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_3 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_4 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_5 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_6 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_7 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_8 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_9 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_10 |	NADH-DEHYDROG-A-RXN |
| S_NS1_Af_002_v1v9_11 |	RXN-11135 |
| S_NS1_Af_002_v1v9_12 |	RXN-11135 |
| S_NS1_Af_002_v1v9_13 |	RXN-11135 |
| S_NS1_Af_002_v1v9_14 |	RXN-11135 |
| S_NS1_Af_002_v1v9_15 |	RXN-11135 |
| S_NS1_Af_002_v1v9_16 |	RXN-11135 |
| S_NS1_Af_002_v1v9_17 |	3.2.1.23-RXN |
| S_NS1_Af_002_v1v9_17 |	RXN-12400 |
| S_NS1_Af_002_v1v9_17 |	RXN-12399 |
| S_NS1_Af_002_v1v9_17 |	RXN-12398 |
| S_NS1_Af_002_v1v9_17 |	BETAGALACTOSID-RXN |
| S_NS1_Af_002_v1v9_17 |	KETOLACTOSE-RXN |

In this example, that genome has10 copies of EC 1.6.5.3, 6 of 3.6.4.12, and one for 3.2.1.23, but the first two ECs are only associated with one reaction each, whereas the last one is associated with 6 reactions. Therefore, the first 16 rows of the output file will have the same reaction ID, and the last row will have 6 different reactions. The report file will contain the ECs that were not present in the PICRUSt2 database, and it will also print the total number of ECs not present for each input file.

Finally, I had to modify the existing PICRUSt2 file where the pathways were listed with the set of reactions but modified by hierarchy. That way the new PICRUSt2's metacyc file is compatible with my version of MinPath. After that, I ran a test with a sample and it lookes like it worked. I will check on it the next working day.

This was the command I ran in the MinPath folder:
`python MinPath.py -any ../MetaCyc/EC_annotated_genomes/Reactions/Sanger/S_NS1_Af_002_minpath_reactions.tsv -map ../MetaCyc/Master_Files/picrust2_pathways_reactions_subset.txt -report any_report -details any_details`


### February 17th, 2026

After generating the new database and new input files for MinPath, I am ready to run my annotated EC genomes with the PICRUSt2's database. The command I ran is this one:
`for file in *; do python ../../../../MinPath_2026/MinPath.py -any ${file} -map /Users/danielcm/Desktop/SickKids/MetaCyc/Master_Files/picrust2_pathways_reactions_subset.txt -report /Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/reports/updated/picrust2_db/${file%%.tsv}_picrust2_report.txt -details /Users/danielcm/Desktop/SickKids/MetaCyc/MetaCyc_Minpath_output/annotated_pathways_genomes/detailed_reports/updated/picrust2_db/${file%%.tsv}_picrust2_detailed_report.txt; done;` 

In this directory: `/Users/danielcm/Desktop/SickKids/MetaCyc/EC_annotated_genomes/Reactions/Barrnap`
And then repeated in this one: `/Users/danielcm/Desktop/SickKids/MetaCyc/EC_annotated_genomes/Reactions/Sanger`

I confirmed that the output files are generated in the output folders from my command. 

After that, I modified and ran the `parsing_minpath_reports.py` script to parse the output files from MinPath and merge them with the Master MetaCyc file to add more information about each pathway. Each folder has different report files, and some of these include the `default` Minpath results (using the default database) and the `updated` MinPath results (using the PICRUSt2's database or Metacyc 2026).

Tomorrow I will be working on comparing the results between the default and updated MinPath runs, and I will also be looking at the differences in pathway annotations between the two approaches. I will be categorizing each pathway as well and generate figures that will simplify my findings for each consortia, method, and database.

### February 18th, 2026

Today I will be working on the categories of my pathways. 

### February 19th, 2026

I had to generate a new file where I am storing the Metacyc categories according to the specific 'Ontology'. This was a complete mess (not new from Metacyc), so I had to do it manually by looking at the 'Ontology - pathway type' column of the Master MetaCyc file and assigning a category to each pathway. This new file is called `metacyc_category_curation_2026.txt` and it has four columns. Each column has a Level 1 or Level 2 value. Level 1 is the highest hierarchy of a category (i.e. Degradation, Biosynthesis, Detoxification, etc.), and level 2 has a more specific category (i.e. Carbohydrate degradation, Amino acid biosynthesis, etc.). The next step is to assign the proper category to my annotated genomes from the minpath, Metacyc 2026, and PICRUSt2 databases. This new modification will be added to my existing python script `parsing_minpath_reports.py`.


### Working with Ayman
Laura wants to try two primers, dec3 and dec5 on the contrived community from Zymo. It has Enterococcus faecalis, Escherichia coli, Lactobacillus fermentum, Listeria monocytogenes, Pseudomonas aeruginosa, Salmonella enterica, Staphylococcus aureus, and Saccharomyces cerevisiae. She wants to try the published primers on this community. We are expecting them to only bind to E. faecalis. We need to download all the strains of E. faecalis from IMG to compare the allelic heterogeneity.

#### Back to T1D

I carried on with the individual annotation of each genome using three databases. In the end, the script that matters the most is `parsing_minpath_reports.py` because it is the one that adds the categories and summaries to each pathway found in the MinPath output files. It can be explained like this:

```
This script parses and summarizes MinPath pathway annotation reports generated from genome files using two annotation methods (default and updated) and two MetaCyc databases (metacyc_2026 and picrust2_db).

What it does
The script is divided into four functions:

sanity_check() — Validates the curated category file by checking for unexpected duplicates across all ontology levels (Level 1 through Level 2.3).

annotating() — Reads raw MinPath .txt report files, extracts naive and parsimonious (MinPath) pathway predictions, and merges them with a master MetaCyc pathway metadata file. Outputs one parsed .tsv file per genome.

simplifying_categories() — Reads each parsed .tsv file and appends ontology classification columns (Level 1–1.3 and Level 2–2.3) by looking up each pathway's ontology type in a curated category file. Categories are sorted alphabetically with NA values placed last.

counting_classification() — Reads all parsed and categorized files and generates summary statistics per genome, including total pathways predicted, naive counts, MinPath counts, bacterial vs non-bacterial classifications, and counts per ontology category.

Input
Raw MinPath report .txt files (one per genome)
Master MetaCyc pathway metadata file (.tsv)
Curated ontology category file (.txt)
Output
One parsed .tsv file per genome with pathway predictions and metadata
Summary classification count files per method and database combination
```

### February 20th, 2026

I had a meeeting with Laura, and I am working on the next set of results and finalized version of figures for the T1D manuscript. I am first re-doing the PICRUSt2 results with some modified parameters for me to get the contribution of ASV, mouse, per pathway per community. I am running this in ARC in `/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/MASTER/Mice/PICRUSt2.3`. Here is the first command:

`for community in {"ns1","ns6","s2","s5"}; do place_seqs.py -s ../Phyloseq/ps_${community}_asv_sequences.fasta -r ${community}_local_files -o ${community}_output/${community}_placed_seqs.tre -p 10 --intermediate ${community}_output/asv_intermediate_seqs_${community} --verbose; done;`

This command produced the following output text per community (showing just NS1):

```
Raw input sequences ranged in length from 279 to 427

epa-ng --tree ns1_local_files/ns1_local_files.tre --ref-msa ns1_output/asv_intermediate_seqs_ns1/ref_seqs_hmmalign.fasta --query ns1_output/asv_intermediate_seqs_ns1/study_seqs_hmmalign.fasta --chunk-size 5000 -T 10 -m ns1_local_files/ns1_local_files.model -w ns1_output/asv_intermediate_seqs_ns1/epa_out --filter-acc-lwr 0.99 --filter-max 100
INFO Selected: Output dir: ns1_output/asv_intermediate_seqs_ns1/epa_out/
INFO Selected: Query file: ns1_output/asv_intermediate_seqs_ns1/study_seqs_hmmalign.fasta
INFO Selected: Tree file: ns1_local_files/ns1_local_files.tre
INFO Selected: Reference MSA: ns1_output/asv_intermediate_seqs_ns1/ref_seqs_hmmalign.fasta
INFO Selected: Filtering by accumulated threshold: 0.99
INFO Selected: Maximum number of placements per query: 100
INFO Selected: Automatic switching of use of per rate scalers
INFO Selected: Preserving the root of the input tree
INFO Selected: Specified model file: ns1_local_files/ns1_local_files.model
INFO Selected: Reading queries in chunks of: 5000
INFO Selected: Using threads: 10
INFO     ______ ____   ___           _   __ ______
        / ____// __ \ /   |         / | / // ____/
       / __/  / /_/ // /| | ______ /  |/ // / __  
      / /___ / ____// ___ |/_____// /|  // /_/ /  
     /_____//_/    /_/  |_|      /_/ |_/ \____/ (v0.3.8)
INFO Using model parameters:
INFO    Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 0.369665 (user),  weights&rates: (0.25,0.0125923) (0.25,0.158759) (0.25,0.696251) (0.25,3.1324) 
        Base frequencies (user): 0.233019 0.241872 0.327604 0.197505 
        Substitution rates (user): 0.758281 1.92043 1.4526 0.755943 3.18076 1
INFO Output file: ns1_output/asv_intermediate_seqs_ns1/epa_out/epa_result.jplace
INFO 179 Sequences done!
INFO Time spent placing: 0s
INFO Elapsed Time: 0s


gappa examine graft --jplace-path ns1_output/asv_intermediate_seqs_ns1/epa_out/epa_result_parsed.jplace --fully-resolve --out-dir ns1_output/asv_intermediate_seqs_ns1/epa_out
                                              ....      ....  
                                             '' '||.   .||'   
                                                  ||  ||      
                                                  '|.|'       
     ...'   ....   ... ...  ... ...   ....        .|'|.       
    |  ||  '' .||   ||'  ||  ||'  || '' .||      .|'  ||      
     |''   .|' ||   ||    |  ||    | .|' ||     .|'|.  ||     
    '....  '|..'|'. ||...'   ||...'  '|..'|.    '||'    ||:.  
    '....'          ||       ||                               
                   ''''     ''''   v0.8.5 (c) 2017-2024
                                   by Lucas Czech and Pierre Barbera

Invocation:                        gappa examine graft --jplace-path ns1_output/asv_intermediate_seqs_ns1/epa_out/epa_result_parsed.jplace --fully-resolve --out-dir ns1_output/asv_intermediate_seqs_ns1/epa_out
Command:                           gappa examine graft

Input:
  --jplace-path                    ns1_output/asv_intermediate_seqs_ns1/epa_out/epa_result_parsed.jplace

Settings:
  --fully-resolve                  true
  --name-prefix

Output:
  --out-dir                        ns1_output/asv_intermediate_seqs_ns1/epa_out
  --file-prefix                    
  --file-suffix

Newick Tree Output:
  --newick-tree-quote-invalid-chars
                                   false

Global Options:
  --allow-file-overwriting         false
  --verbose                        false
  --threads                        16
  --log-file

Run the following command to get the references that need to be cited:
`gappa tools citation Czech2020-genesis-and-gappa`

Started 2026-02-20 16:19:23

Found 1 jplace file

Finished 2026-02-20 16:19:23
```

The next part was to run the hidden state prediction placement. For that I ran the following command:
`for community in {"ns1","ns6","s2","s5"}; do hsp.py -t ${community}_output/${community}_placed_seqs.tre --observed_trait_table trait_tables/16S.txt --calculate_NSTI -p 10 --seed 23 --verbose -o ${community}_output/${community}_16S_nsti.predicted.tsv; done`

which gave me the following output:

```
Rscript /home/daniel.castanedamogo/anaconda3/envs/picrust2/lib/python3.9/site-packages/picrust2/Rscripts/castor_nsti.R ns1_output/ns1_placed_seqs.tre /home/daniel.castanedamogo/tmp/tmpsold53vi/known_tips.txt /home/daniel.castanedamogo/tmp/tmpsold53vi/nsti_out.txt

Warning messages:
1: package ‘castor’ was built under R version 4.4.1 
2: package ‘Rcpp’ was built under R version 4.4.1 
```

I repeated the same command across the KOs and ECs trait tables, and I got the predicted tables for each community using this command:

`for community in {"ns1","ns6","s2","s5"}; do for function in {"KO","EC"}; do hsp.py -t ${community}_output/${community}_placed_seqs.tre --observed_trait_table trait_tables/${function}_for_picrust2_renamed.tsv --calculate_NSTI -p 10 --seed 23 --verbose -o ${community}_output/${community}_${function}_nsti.predicted.tsv; done; done;`

 The next step is to run the metagenome prediction:

`for community in {ns1,ns6,s2,s5}; do for function in {KO,EC}; do metagenome_pipeline.py -i ${community}_input/ps_${community}_final.biom -m ${community}_output/${community}_16S_nsti.predicted.tsv -f ${community}_output/${community}_${function}_nsti.predicted.tsv --strat_out --wide_table -o ${community}_output/${community}_${function}_metagenome_out; done; done;`

This time, I am including the `--wide_out` and `--strat_out` parameters to get the output in a wide format and to get the stratified output as well. The wide output will be easier to work with for me, and the stratified output will allow me to see the contribution of each ASV to each pathway.

and I got all the results under a 2.0 of NSTI, which is the recommended threshold for PICRUSt2. This is expected, the lowest NSTI is tiny.

```
All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.
All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.
```

The stratified output looks like this:
```
function	sequence	Plate1_1030R_0_M_NS1_week10_S77_L001
EC:1.6.5.3	ASVp1_1	97295.38
EC:1.6.5.3	ASVp1_10	3648.75
EC:1.6.5.3	ASVp1_104	120.60000000000001
EC:1.6.5.3	ASVp1_12	3196.0
```

Where I have the normalized abundance and contribution of each ASV to each pathway. The unstratified output file does not provide contribution of each ASV, and it only goes from EC to EC.

Finally, I am running the `pathway_pipeline.py` to get the pathway abundance predictions from the ECs. This is the command I am running:

`for community in {ns1,ns6,s2,s5}; do pathway_pipeline.py -i ${community}_output/${community}_EC_metagenome_out/pred_metagenome_strat.tsv -o ${community}_output/${community}_pathway_inference/ --intermediate ${community}_output/${community}_pathway_intermediate_files/ -p 32 --coverage --per_sequence_contrib --per_sequence_abun ${community}_output/${community}_EC_metagenome_out/seqtab_norm.tsv --per_sequence_function ${community}_output/${community}_EC_nsti.predicted.tsv --wide_table --verbose; done;`

It produces a bunch of output files:

- `path_abun_strat.tsv.gz`
   Option to specify that stratified abundances should be reported in terms of the contribution by each predicted genome rather than how much each genome is contributing to the overall community. In other words, pathway abundances will be calculated for each individual predicted genome. Both --per_sequence_abun and --per_sequence_function need to be specified when this option is set. Stratified coverages will only be reported when this option is used (and --coverage is set). As of v2.2.0-b, unstratified pathway abundances based on the community-wide pathway abundances and also based on the per-sequence pathway abundances will be output when this option is used.

- `path_abun_predictions.tsv`
ASV contribution to the overall community pathway abundance (community, not mouse or sample). 

- `path_abun_strat.tsv`
Pathway abundance per ASV per sample.

- `path_abun_unstrat.tsv`
Pathway abundance per sample, all ASVs contributions per sample are summed, no ASV breakdown.

- `path_abun_unstrat_per_seq.tsv`

- `path_abun_unstrat_per_seq.tsv`
ASV contribution to total pathway abundance across samples. This file is different from path_abun_unstrat.tsv in how the pathway abundances were calculated. This filecomputes pathway abundance by considering individual ASV contribution, whereas the 'path_abun_unstrat.tsv' only collapses to ECs first without considering ASVs at all.

`path_cov_predictions.tsv` 
Predicts if a pathway is encoded per ASV. Each number here is between 0 - 1. That score assesses how likely a pathway is present, the closer to 1, the more likely it is produced by that ASV.

- `path_cov_strat.tsv` Same idea, but in here it is per ASV per sample.

- `path_cov_unstrat.tsv` Same idea, but in here it is only pathway core per sample.

### February 21st, 2026

Today I will use the picrust2 files as input and crossmap each ASV to their taxa using the Phyloseq file I have. For that, I created the script `asv_to_taxa_for_picrust2.3.py` that takes the ASV IDs from the PICRUSt2 output files and matches them to the corresponding taxa in the Phyloseq files. I have uploaded the new files into the PICRUSt2 folder of my SickKids repository and the One Drive too.

I managed to make another script that makes a boxplot figure of the averages of a level 2 categories that are enriched in maaslin2 for either side of the Log2(FC). This script is called `enriched_maaslin2_category_average.R`. I will make a new script called `bubble_plot_taxa_pathways_mice.R` that will attempt to plot a bubble plot with many variables. 

### February 26th, 2026

We realized there's a bug in 'Have-some-pie', where the merging of the ASVs start pointing at different taxa. I will be fixing that today by looking at `phyloseq_t1d_db2.1`. First, I need to remember that 4 samples are not labeled with their week, and those are:

```
Plate4_036R_0_M_NS1_S62_L001,
Plate4_036R_L_M_NS1_S63_L001,
Plate4_036R_RL_M_NS1_S64_L001,
Plate4_036R_RR_M_NS1_S65_L001,
```

which are found in `/Users/danielcm/Desktop/diammatics/T1D`. I will be re-labeling these four samples in the plate 4 of the `/Users/danielcm/Desktop/Sycuro/Projects/Diabetes/t1d_db_fixed_discussed/FemMicro_Daniel/` folders for plate1.1, plate2.1, plate3.1, plate4.1, and plate5.1. 

It looks like I found the culprit in my phyloseq_code, which was in the 'writing_ps_objects' (or something along those lines). I removed it and simplify my code even further. Its name will remain the same: `phyloseq_t1d_db2_2.1.R`. I will be running this code across all the plates to generate new phyloseq objects with the correct sample names. 

I was able to run the first two steps of PICRUSt2: place_seqs and hsp.py. Tomorrow I will resume with the metagenome_pipeline.py.


### March 1st 2026

Laura found an error in taxonomy assignment in the latest files she had access to. Luckily, we realized that our collaborators had the right files all along, so it didn't mess up any important analysis. Regardless of that, I decided to find out why were we short by ~550 ASV counts out of 23 Million. Turns out I used the wrong non_chimeric_object.rds from plate1 instead of plate1.1. I decided to run the rest of PICRUSt2, Maaslin2, and others just to make sure we have the proper results.

### March 6th, 2026

I had a discussion with the team after we had to review a table that 'had contaminants'. Regardless of those contaminants not being present at the prevalence that were shown to us, a filtering step was requested to remove ASVs that are in low abundance and are contaminants. I will be doing this on the mice table and the consortium table. I will be calculating the sensitivity (TP/TP+FN) and the specificity (TN/TN+FP).

### March 9th, 2026

Today I will focus on the filtering approaches across the inocula samples as my gold standard, and the mice samples after finding the best threshold for filtering.

NS1:
Using the consortia NS1 as our reference, I realized that we need a minimum of 173 ASV counts to remove the contaminant and without loosing any true positives. It's equivalent proportion is 0.1237%, and it is present in all 3 samples. If I apply a filter of 175 ASV count, then it is removed, and the rest of the true taxa are kept. False positive 1 to 0

S2:
There are three contaminants or rare ASVs (Enterococcus faecalis, NA NA, NA NA), where the max sum for each is 3,4, and 11. By simply applying the previous filter of 175, we remove those contaminants, but we also remove a true positive for the only ASV pointing at Collinsella sp902362275. Collinsella is present in all 3 samples. False negative = 1, false positive 3 to 0

NS6:
In NS6, there are plenty ASVs pointing at taxa that are not supposed to be part of that consortium. Most of them are linked to the new Cytobacillus bacterium. The highest ASV count for a contaminant is 653 for A. muciniphila, followed by 231 for B. thetaiotaomicron, and the 3rd one is a contaminant from Flavobacterium ammonificans with 164. If we apply the 175 ASV count rule, we fail to remove A. muciniphila, and B. thetaiotaomicron, but we remove A. faecis, CHH4-2, Cloacibacterium sp902362275, Cytobacillus spp. E. alcoholdehydrogenati, adding up to 14 taxa that are removed. 

False positives from 18 to 4
False negatives from 0 to 2

S5:
When applying the 175 ASV filter, we remove 4 false positives (B. stercoris, B. uniformis, P. merdae, S. wadsworthensis), we do not lose any true positives this case)


PICRUSt2 workflow:
1. Gene content inference is calculated for each organism based on the phylogeny tree.
2. PICRUSt2 may predict gene presence on organisms that have not been sequenced yet based in their sequenced evolutionary relatives.
3. PICRUSt2 requires the 16S copy number for each organism being analyzed, and to do that, it looks up its value in the reference database.
4. The user's ASV abundance table is normalized by its predicted 16S copy number.
5. Normalized ASV abundances are then multiplied by the set of gene family abundances pre-calculated in the hidden state prediction step (ancestral state reconstruction) to yield the final metagenome prediction.

### March 12th, 2026

Now I will apply the 175 ASV rule to the rest of the mice groups to see how 'good' this is for removing FPs while keeping TPs.

As a reminder for later: NS1 mice = 93, NS6 mice = 32, S2 mice = 106, S5 mice = 32.


### March 13th, 2026
I've been tunning new parameters to minimize the number of false positives as much as I can. I have generated a new code called `phyloseq_2.2.R` that tunes and iterates parameters to find the best filter values across each consortia. I will do a sanity check and debug if needed, but so far it looks promising.

### March 16th, 2026
I have generated a new code called `phyloseq_2.2.R` that tunes and iterates parameters to find the best filter values across each consortia. I will do a sanity check and debug if needed, but so far it looks promising. The results I have so far are the following:

Best parameters for NS1 (FP = 0, Not losing any ASVs pointing at the true positives):
- ASV count of at least 600 or more
- ASV length of 200 or more (can be extended to 400)
- ASV prevalence across samples of 0.10% or more (can be extended to 1%)

with 1 FP, this can be modified to:
- ASV count: 200 or more (can be extendedto 550 bp)
- ASV length: 200 or more (can be extended to 400)
- ASV prevalence across samples of 0.1% or more (can be extended to 1%)



Best parameters for NS6 (FP = 0, Not losing any taxon pointing at the true positives):
- ASV count of at least 200 or more (can be extended to 300)
- ASV length of 200 or more (can be extended to 400)
- ASV prevalence across samples of 5% or more

With 1 FP, this can be modified to:
- ASV count: 100 or more (can be extended to 500)
- ASV length: 200 or more (can be extended to 400)
- ASV prevalence across samples of 0.1% (can be extended to 1%)



Best parameters for S2 (FP = 0, not losing any taxon) 
- ASV count of at least 300 or more (can be extended to 600)
- ASV length of 200 or more (can be extended to 400)
- ASV prevalence of 15% or more (can be extended to 20%)

With 2 FP, this can be modified to (it jumps from 0 to 2):
- ASV count of at least 350 or more (can be extended to 600)
- ASV length of 200 or more (can be extended to 400)
- ASV prevalence of at least 10% or more (cannot be changed)


Best parameters for S5 (FP = 2, this is the lowest):
- ASV count of at least 300 or more (can be extended to 600)
- ASV length of at least 200 or more (can be extended to 400)
- ASV prevalence of at least 20% or more

Because its FP is already too high, I don't recommend changing these parameters at all.

Now I will test the hsp.py prediction from Castor by taking some ASVs that seem relatively distant from the rest. For instance, taking this one from the engrafted tree produced by epa-ng and gappa:

or taking this one `TGGGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGCAGGATGAAGGTCTTCGGATTGTAAACTGCTTTTGTCAGGGACGAAAAGGGATGCGATAACACCGTATTCCGCTGACGGTACCTGAAGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTCTGTAAGATAGATGTGAAATCCCCGGGCTCAACCTGGGAATTGCATATATGACTGCAGGACTTGAGTTTGTCAGAGGAGGGTGGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAAGAACACCGATGGCGAAGGCAGCCCTCTGGGACATGACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACA`

both FemMicro and BLAST agree that this is S. wadsworthensis, which is supposed to be in NS6 (and in NS1 too). In the tree (see below if I succeeded at copying the screenshot), this ASV is closest to Sw 061, and its next 'known' parent would be Ef 

When looking at the KO and EC tree from this ASV vs from the annotation profile I did on the genome for S. wadsworthensis, I get these metrics:

Total ECs across all genomes: 1,886
ECs with no difference in count between ASV and E.faecalis: 1,230 (65.2%)
ECs with no difference in count between ASV and S. wadsworthensis: 1,665 (88.3%)

This decides it. It is clear that the profile is not identical and even after playing with penalization numbers the annotation numbers do not identically overlap. I will be working on bypassing this step from PICRUSt2 and just using the ECs that I have annotated for each genome. This way, I can be sure that the profiles are identical and that the differences in pathway predictions are due to the differences in the databases and not due to the hidden state prediction step from Castor, the name of the script will be `hsp_modified_castor.py`


### March 17th, 2026

I will be including now the False Negative data as part of the analysis. And prepare my data update for tomorrow, my last one, probably.

### March 23rd, 2026

Last week we had a meeting, and in there we evaluated the results of the filtering step, and we decided to apply the following filters to the mice and inocula samples, across each inoculum group, keeping the ASVs that fall in these criteria:

- ASV length of 400 bp or more (inclusive) (400 bp was recently changed from 250 bp)
- ASV total count of 600 or more (inclusive)
- ASV prevalence across samples of 0.20 (20%) or more (inclusive)

I tested these new parameters to see how the EC count profile difer from the genomes (the script I made for this is called `castor_testing.py`), and I noticed that the biggest factor decreasing dissimilarities in EC profile is not the edge_cost of castor itself (set to 0, 0.5, 2, 5, and 10), but applying the ASV filter. When setting the penalty to zero in the `hsp.py` script from PICRUSt2, we get a higher similarity in EC profile count bewteen ASVs and ECs from the genome. This is expected, as the cost formula between EC transitions is: `cost = 1/(edge_length^edge_cost)`, so when edge_cost is set to zero, the cost of transition between ECs is 1, which means we are ignoring the length of the branch in the tree and going for the nearest neighbour common ancestor.

I took the 220 16S sequences (barrnap + sanger) that have some duplicates, and using `seqkit rmdup`, I figured out which sequences are duplicated. After that, I modified these repeated sequences by modifying the first or 2nd nucleotide in the sequence. Then I ran a `mafft` alignment on those, `hmmbuild` and put them in the `all_consortia_reference_files` folder, located in the `MASTER/PICRUSt2.4` directory.


### March 30th, 2026

I tried a new approach by modifying the code from PICRUSt2, so it would accept Castor's method of Nearest Neighbour Joining. In a nutshell, this is what I did:

1. Created a new conda environment by clonning my existing one: `conda create --name picrust2_mod --clone picrust2.2`
2. I activated the new environment: `conda activate picrust2_mod`
3. I wanted to see where the `hsp.py` script is located, so I ran `which hsp.py`, and it showed me the path: `~/anaconda3/envs/picrust2_mod/bin/hsp.py`, this is telling me that this code is in the bin folder of the conda environment, but it should be called by a wrap script
4. I wanted to be extra safe and I made a PICRUSt2 backup: `conda env export > picrust2_backup.yml`
5. I wanted to see where the picrust2 scripts are located, so I ran `python -c "import picrust2, os; print(os.path.dirname(picrust2.__file__))"`, and it showed me the path: `/home/daniel.castanedamogo/anaconda3/envs/picrust2_mod/lib/python3.6/site-packages/picrust2` so now I don't need to be looking for the picrust2 files all over my Anaconda folder
6. I also made new backups of the scripts I am planning on modifying: 
   `cd /home/daniel.castanedamogo/anaconda3/envs/picrust2_mod/lib/python3.6/site-packages/picrust2`
   `cp wrap_hsp.py wrap_hsp.py.bak`
   `cp Rscripts/castor_hsp.R Rscripts/castor_hsp.R.bak`
7. I modified the `Rscripts/castor_hsp.R` to include a new alternative that has the nearest neighbour, so after this part from the original script:
```
 else if (hsp_method == "subtree_average") {

    predict_out <- lapply(trait_values,
                            hsp_subtree_averaging,
                            tree = full_tree,
                            check_input = check_input_set) }
```

I added this:

```
else if (hsp_method == "nearest_neighbor") {

    predict_out <- lapply(trait_values,
                            hsp_nearest_neighbor,
                            tree = full_tree,
                            check_input = check_input_set)
  }
```

I knew I had to add it to this particular section of the script because NN returns predicted states directly, and not likelihoods like `mp` does or `emp-prob` (which have a few variables after calling it). It made sense to put it into the if-else block of `subtree-average` that gives the states directly.

8. I also needed to modify that same script so if the user adds the argument 'nearest_neighbor' it would be recognized. So I looked for the args part of the script and I modified it like this:

```if (hsp_method == "pic" || hsp_method == "scp" || hsp_method == "subtree_average" || hsp_method == "nearest_neighbor") {
   .
   .
   . }
```
9. In the `bin/hsp.py` script, I modified the `HSP_METHODS` so it can include the one I just added, along with a description of what it does:


```
HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average', 'nearest_neighbor']
```

and this:

```
parser.add_argument('-m', '--hsp_method', default='mp',
                    choices=HSP_METHODS,
                    help='HSP method to use.' +
                    '"mp": predict discrete traits using max parsimony. '
                    '"emp_prob": predict discrete traits based on empirical '
                    'state probabilities across tips. "subtree_average": '
                    'predict continuous traits using subtree averaging. '
                    '"pic": predict continuous traits with phylogentic '
                    'independent contrast. "scp": reconstruct continuous '
                    'traits using squared-change parsimony (default: '
                    '%(default)s). '
                    '"nearest_neighbor": method introduced by Daniel CM for '
                    'discrete trait prediction by using nearest neighbor')
```
10. Finally, I ran the new modified scripts with the new method, and it worked! I got the predicted EC table for each ASV. I took that table to my local computer to run it with my script `castor_testing.py` to see how similar the predicted EC profile, and turns out that the EC count profile for the filtered `ns6` table had 0 differences from the genome we expect it is coming from! This worked! Now I need to repeat the same across the rest of the communities, and do it with the unfiltered data and filtered as a sanity check.              


### April 1st, 2026

I ran the picrust2 pathway prediction, and it seems that the ASVs pointing at the same taxa have the same abundance for each pathway, which is expected because they have the same EC profile. I ran this test by looking at Alessandra's value with discrepant data for P. distasonis:

![Alt text](/Users/danielcm/Desktop/ASV_taxa_pathway_test.png)

The next part was to run `merge_functional_predictions.py` to merge the pathway predictions across the filtered resullts, so I have only one file with all of the consortia and all the pathways.

The next step is to run `maaslin2_t1d_db2.R` to generate the enrichment prediction across groups, I had to modify the metadata file so the sample names match the ones Laura coded for. I also had to modify the output file name, and the input directory where the predicted pathway merged file is located.

### April 13th, 2026 (back from vacation)

It turns out some of the IDs that were generated by Stata (...) had some errors, so I need to run them again with the correct IDs. For this, I will flag the ones that are not supposed to be in each consortium and remove the w29 that was grouped into week 9. After cleaning up the ASV tables, I will make a script that changes those names to the ones submitted to NCBI. I will run both tables with both nomenclatures in case we need one over the other.

In parallel, I will be running the Parkinson's Data for my talk next friday (24th, April). 

I have finished up the script `renaming_sample_IDs.py` that takes the old sample names and replaces them with the new ones. I will be running this script across all consortia so it gets the new SRA names. This is now in the repo. The files I have generated with the new names are in `/Phyloseq2/Filtered_to_use/` and they have the SRA ID in the name of the file. These are good ASV count tables and are ready to be analyzed by PICRUSt2 after I create new fasta files for each ID.

### April 14th, 2026

In PICRUSt2.6 I will have the newest results with the fixed IDs. This is what I have done so far:

1. Run `Phyloseq_2.2.R` to generate new .csv files with the filter criteria we decided on. Then I manually transformed the output file of that script (i.e. `NS1_filtered_ASVs_count600_len400_prev20.csv`) to only keep the asv_id column as well as the sample names (one per column). This file has the same nomenclature but ends with `_f.csv` to indicate it is the filtered version. I did this for all the consortia, and I have them in the `Phyloseq2/Filtered_to_use/` folder. 
2. I ran the `renaming_sample_IDs.py` script to change the sample names in the filtered files to the SRA IDs. I have these new files in the same folder, but they end with `_f_sra.csv` to indicate they are filtered and have the SRA IDs.
3. I generated the new fasta files for each consortium with the new ASV IDs, and I have them in the same folder, but they end with `_f_sra.fasta` to indicate they are filtered and have the SRA IDs. To make sure that the names matched between these sample IDs and the metadata, I simply took the IDs of my file and the metadata and put them into Venny, making sure the intersection of both was 100%
4. I copied these files into ARC in the `PICRUSt2.6` folder, and using a bash command, I transformed the `_f_sra.csv` files into `.tsv` by doing: `tr ',' '\t' < NS1_filtered_ASVs_count600_len400_prev20_f_sra.csv > NS1_filtered_ASVs_count600_len400_prev20_f_sra.tsv` and I did this for all the consortia. If I don't, then PICRUSt2 will not recognize the .csv files as an ASV input table.
5. `conda activate picrust2_mod` to activate the modified conda environment with the new hsp method.
6. I ran the `place_seqs` step of PICRUSt2 for each consortium:


