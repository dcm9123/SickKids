# T1D Project Entries from 2026

### February 2nd, 2026

I have decided to create a new repository that only contains the publicly available files as well as the necessary ones (as long as they are below the capacity of Github's storage). I was asked to do an annotation profile of each genome with the ECs and therefore the MetaCyc pathways. After dealing so long with MetaCyc, Kevin and I figured out that the best way to do it is to export the 'smart tables' by indicating which fields we are interested in. However, the website is crap and we can only do 5 fields at a time, then I came up with a script `pathway_metacyc_merger.py` that takes the fields I need across multiple files, and then they are merged into one large one, called `Master_Metacyc_pathway_file.tsv`. I had to come up with a list of bacterial names, genus, family, etc in order to keep a database exclusively for bacteria, as their fields are really bad and not filtered by Kingdom.

I have finished coding the `metacyc_genome_annotator.py` which retrieves the genome number and prints the EC multiple times (n times of CNV found for that cell). The input file is the `EC_for_picrust2_renamted.tsv` and it generates an output file for every row, printing the number and ID of each EC. I made sure all rows had a different ID by simply adding a counter. Now that I have both files, MinPath should run smoothly. But that's a task for tomorrow!

### February 3rd, 2026

I ran the `MinPath` software across one sample to test it out. This is the command I used: `python MinPath.py -ec ../../diammatics/T1D/PICRUSt2.2/EC_annotated_genomes/S_NS1_Af_002_minpath_ecs.tsv -map ../Minpath_ready_Metacyc_pathway_file.tsv -report test_report -details test_details` inside this directory: `/Users/danielcm/Desktop/SickKids/MinPath_2026`

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

Here is the difference between the naïve and the parsimonious approach in a figure *that does not belong to me*:
https://www.google.com/url?sa=t&source=web&rct=j&url=https%3A%2F%2Fjournals.plos.org%2Fploscompbiol%2Farticle%3Fid%3D10.1371%2Fjournal.pcbi.1000465&ved=0CBYQjRxqFwoTCPDrnfSyvpIDFQAAAAAdAAAAABBO&opi=89978449<img width="600" height="407" alt="image" src="https://github.com/user-attachments/assets/433f6a0f-9d34-4e39-9573-a0ccbfc8ab38" />


