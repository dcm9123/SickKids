# T1D Project ENTRIES FROM 2026 ONLY

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

