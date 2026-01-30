# Daniel Castaneda Mogollon, PhD
# January 29th, 2026
# Script designed to merge multiple MetaCyc pathway, reaction, and enzyme files into a single file for easier processing

# %%
import pandas as pd
import os
import re

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# %%
def filtering_by_bacteria():
    
    os.chdir("/Users/danielcm/Desktop/diammatics/T1D/MetaCyc_files/")
    df = pd.read_csv("pathways_with_species.txt", sep="\t")
    df['Common-Name'] = df['Common-Name'].str.strip()
    df['Common-Name2'] = df['Common-Name2'].str.strip()

    df['Common-Name'] = df['Common-Name'].str.replace(' // ',',')
    df['Common-Name2'] = df['Common-Name2'].str.replace(' // ',',')
    
    
    regex_bacteria = re.compile(
    r"""(?ix)                           # i=ignorecase, x=verbose
    (?:\bBacteria\b|<bacteria>)         # explicit bacteria markers

    | \b(?:                             # higher taxonomic bacterial groups seen in your file
        Pseudomonadota|Actinomycetota|Bacillota|Cyanobacteriota|
        Gammaproteobacteria|Alphaproteobacteria|Betaproteobacteria|
        Bacteroidales|Lactobacillales|Campylobacterales|
        Enterobacteriaceae|Prevotellaceae|Microbacteriaceae|Chlorobiaceae|Alcaligenaceae|
        Streptomycetaceae|Micromonosporales|Kitasatosporales|
        Clostridia
    )\b

    | \b(?:                             # bacterial genera (expanded to include ones you actually have)
        Escherichia|Salmonella|Pseudomonas|Streptomyces|Bacillus|Clostridium|Clostridioides|
        Enterococcus|Klebsiella|Acinetobacter|Rhodococcus|Cupriavidus|Mycobacterium|Corynebacterium|
        Vibrio|Yersinia|Azotobacter|Bifidobacterium|Lacticaseibacillus|Lactiplantibacillus|Phocaeicola|
        Bacteroides|Prevotella|Paenibacillus|Campylobacter|Helicobacter|Bordetella|Achromobacter|
        Rhodobacter|Chlorobaculum|Chlorobium|Allochromatium|Methylobacillus|
        Agrobacterium|Microbacterium|Hoylesella|Xylanibacter|Adlercreutzia |Neisseria|Shewanella|Nitrosomonas|
        Haemophilus|Brucella|Rhizobium|Shinorhizobium|Chlamydia|Aquifex|Hydrogenobacter|Prochlorococcus|Chloroflexux|Myxococcus|Stigmatella|
        Shingomonas|Sphingobium|Novosphingobium|Candydatus|Burkholderia|Paraburkholderia|Chloroflexus|Roseiflexus|Blasthochloris|Methylococcus|
        Streptococcus|Staphylococcus|Listeria|Enterobacter|Citrobacter|Serratia|Providencia|Morganella|Proteus|Xanthomonas|Ralstonia|
        Mariprofundus|Mycoplasmoides|Desulfolutivibrio|Vibrio|Nitratidesulfovibrio|Oleidesulfovibrio|Fischerella|Kitasatospora|Thiocapsa|
        Porphyromonas|Tannerella|Parabacteroides|Alistipes|Eubacterium|Collinsella|Eggerthella|Anaerostipes|Dialister|Candidatus|
        Marichromatium|Micromonospora|Blastochloris|Chloroflexus|Corallococcus|Desulfocapsa|Streptantibioticus|Treponema|Plesiocystis|
        Enhygromyxa|Methylorubrum|Aromatoleum|Thauera|Ruegeria|Roseovarius|Cereibacter|Fulvimarina|Maritimibacter|Pseudooceanicola|Roseibium|Yoonia|Sulfitobacter|
        Oceanimonas|Halodesulfovibrio|Alcaligenes|Methylotuvimicrobium|Rubrivivax|Bradyrhizobium|Geobacter|Xanthobacter|Oleomonas|Bilophila|
        Synechocystis|Synechococcus|Gloeobacter|Picosynechococcus|Caldanaerobacter|Pseudoramibacter|Lactococcus|Peptoniphilus|Propionibacterium|
        Polaribacter|Desulfovibrio|Pseudoalteromonas|Zymomonas|Thermotoga|Pelobacter|Acetobacterium|Acidiphilum|Acidithiobacillus|Advenella|Aeropyrum|
        Aliivibrio|Alkalicoccus|Alycycliphilus|Barkera|Castellaniella|Cellulomonas|Comamonas|Curtobacterium|Caldicellulosiruptor|Dechloromonas|
        Halomonas|Pantoea|Acidovorax|Akkermansia|Afipia|Algoriphagus|Aminobacter|Arthrobacter|Brevibacterium|Chromobacterium|Cutibacterium|Dickeya|Erwinia|Ectopseudomonas|Erythrobacter|
        Gemmata|Gimesia|Gluconobacter
    )\s+[a-z][a-z-]+(?:\s+\w+)?\b       # genus species (optional strain/extra token)
    """,
    re.IGNORECASE)
    
    regex_non_bacteria = re.compile(
    r"""(?ix)
    \b(?:                               # broad non-bacteria domains/kingdoms/clades
        Archaea|Eukaryota|Viridiplantae|Embryophyta|Magnoliopsida|
        Metazoa|Mammalia|Chordata|Vertebrata|
        Fungi|Ascomycota|Basidiomycota|Opisthokonta|
        Chlorophyta|Bryophyta|Stramenopiles|Alveolata|Amoebozoa|
        Apicomplexa|Euglenozoa|Haptophyta|Bacillariophyta|Oomycota|
        Viruses|Virus|Viroid
    )\b

    | \b(?:                               # archaeal groups 
        Methanobacteria|Methanomicrobia|Methanopyri|Archaeoglobi|
        Methanobacteriales|Methanomicrobiales|Methanosarcinales|Haloferax|Haloarcula|Halobacterium|
        Halalkalicoccus|Halomicrobium|Haloterrigena|Halorubrum|Natrialba|Ignococcus|Natromonas|Metallosphaera|Pyrobaculum|
        Pyrococcus|Sulfolobus|Saccharolobus|Thermoproteus|Methanolobus|Hyperthermus|Aspergillus|Penicillium|Fusarium|Trichoderma|Botrytis|
        Claviceps|Epichloe|Aureobasidium|Armillaria|Mortierella|Neurospora|Schizosaccaromyces|Kluyveromyces|Rhodotorula|Moniliella|Starmerella|
        Yarrowia|Ogataea|Komagataella|Pestalotiopsis
    )\b

    | \b(?:                               # archaeal/fungal/plant genera
        Methanobrevibacter|Methanosarcina|Methanocaldococcus|Methanococcus|Methanococcoides|
        Methanothermobacter|Archaeoglobus|Pyrococcus|Pyrobaculum|
        Sulfolobus|Saccharolobus|Metallosphaera|Thermoproteus|
        Haloferax|Haloarcula|Halobacterium|Halomicrobium|Haloterrigena|Halorubrum|
        Aspergillus|Penicillium|Fusarium|Saccharomyces|Candida|Neurospora|
        Homo|Mus|Arabidopsis|Thermococcus|Thermoplasma|Desulfurococcus|Staphylothermus|Acidianus|
        Amanita|Beauveria|Botryococcus|Chaetomium|Cladonia|Hapsidospora|Mycena|Neonothopanus|Omphalotus|
        Paxillus|Suillus|Tapinella|Tolypocladium|Talaromyces|Mycosarcoma|Phytophtora|Globisporangium|Pseudohyphozyma|
        Camellia|Coffea|Glycine|Oryza|Triticum|Vitis|Solanum|Brassica|Pisum|Medicago|Lupinus|Trifolium|Helianthus|Gossypium|
        Zea|Cucumis|Carthamus|Spinacia|Apium|Phaseolus|Rosa|Perilla|Lavandula|Dendrobium|Cymbidium|Catharanthus|Petroselinum|
        Foeniculum|Capsicum|Lactuca|Malus|Prunus|Fragaria|Citrus|Nicotiana|Populus|Salix|Eucalyptus|Quercus|Fagus|Betula|Juglans|
        Chrysosplenium|Cicer|Coptis|Cerastium|Stellaria|rynchos|Crysanthemum|Dahlia|Pericallis|Lampranthus|evia|Colchium|Cupressus|
        Mentha|Carum|Schizonepeta|Humulus|Pinus|Hypericum|Panax|Litchi|Petunia|Lathyrus|Ricinus|Adonis|Albizia|Ginkgo|lvia|Lithosperum|
        Rubia|Anemone|Curcuma|Zingiber|Adiantum|Isatis|Lepidum|Sinapsis|Hevea|Linum|Lens|Vigna|Leucaena|Narcissus|Galanthus|Lycoris|Nerine|Zephyranthes|
        Impatiens|Lawsonia|Triglochin|Canavalia|Dianthus|Allium|Glycyrrhiza|Piper|Digitalis|Avena|Cinchona|Clarkia|Lotus|Manihot|Carapichea|Alangium|Rosa|
        Hordeum|Veratrum|Dalbergia|Picea|Abis|Amorpha|Podophyllum|Rauvolfia|Cannabis|Aquatica|Limonium|Sophora|Centaurium|Anchusa|Portulaca|Coleus|Tephrosia|Sorghum|
        Crocus|Gardenia|Plectranthus|Jacquinia|Ephedra|Plumbago|Berberis|Eschscholzia|Papaver|Vicia|Davallia|Cydia|Ricordea|Dipsastraea|Renilla|Psilocybe|Stereocaulon|
        Digenea|Palmaria|Fucus|Laminaria|Spodoptera|Watasenia|Strychnos|Stevia|Colchicum|Dictyostelium|Arenicola|Ascaris|Fasciola|Magallana|
        Mytilus|Sipunculus|Salvia|Lithospermum|Apis|Caenorhabditis|Drosophila|Dugesia|Periplaneta|Dactylopius|Aedes|Blaberus|Blatta|Byrsotria|Cyperus|Diploptera|Locusta|
        Nauphoeta|Rhyparobia|Corbicula|Ips|Anemonia|Discosoma|Galleria|Musca|Doryteuthis|Arctia|Diacrisia|Estigmene|Apantesis|Tyria|Pseudo-nitzschia|Chrysopa|
        Myzus|Acheta|Cadra|Lipomyces|Anadara|Arabella|Busycotypus|Cellana|Halichondria|Haliotis|Littorina|Marphysa|Meretrix|Patiria|Pecten|Suberites|
        Zoanthus|Aequorea|Euglena|Plasmodium|Epiactis|Giardia|Santalum|Viburnum|Crithidia|Gallus|Chondrus|Anisodus|Atropa|
        Datura|Hyoscyamus|Alexandrium|Lingulaulax|Noctiluca|Protoperidinium|Pyrocystis|Taxus|Bemisia|Bombyx|Carcinus|Manduca|Penaeus|Faxonius|Schistocerca|Onchidium|Todarodes|
        Glebionis|Rudbeckia|Tagetes|Anethum|Ocimum|Pimpinella|Erythroxylum|Nepeta|Rubus|Fagopyrum|Rhus|Aloe|Tanacetum|Ruta|Petiveria|Conium|Dryopteris|Citrus|Abies|Camptotheca|
        Nothapodytes|Ophiorrhiza|Artemisia|Asparagus|Cryptomeria|Astragalus|Beta|Daucus|Bruguiera|Rhizophora|Kandelia|Vanilla|Lepidium|Rhododendron|Musa|Anigozanthos|Wachendorfia|
        Antirrhinum|Coreopsis|Cistus|Persicaria

        
    )\s+[a-z][a-z-]+(?:\s+\w+)?\b
    """,
    re.IGNORECASE)
    # Add na=False to all .str.contains() calls
    bacteria_present1 = df['Common-Name'].str.contains(regex_bacteria, na=False)
    bacteria_present2 = df['Common-Name2'].str.contains(regex_bacteria, na=False)
    no_bacteria1 = df['Common-Name'].str.contains(regex_non_bacteria, na=False)
    no_bacteria2 = df['Common-Name2'].str.contains(regex_non_bacteria, na=False)
         
    #print(bacteria_present1)
    
    bacteria_dict = {}
    for i,row in df.iterrows():
        pathway = row['Pathways']
        #Case 1, bacteria present in either Common-Name or Common-Name2, regardless of the presence of Eukarya/Archaea
        if (bacteria_present1.iloc[i] == True) or (bacteria_present2.iloc[i] == True):
            bacteria_dict[pathway] = {"Index": i, "Common-Name": row['Common-Name'], "Common-Name2": row['Common-Name2'], "Classification": "Bacteria"}
        #Case 2, Eukarya/Archaea present in either Common-Name or Common-Name2 and Bacteria absent in both
        elif ((no_bacteria1.iloc[i] == True) or (no_bacteria2.iloc[i] == True)) and ((bacteria_present1.iloc[i] == False) and (bacteria_present2.iloc[i] == False)):
            bacteria_dict[pathway] = {"Index": i, "Common-Name": row['Common-Name'], "Common-Name2": row['Common-Name2'], "Classification": "Non-Bacteria"}
        # Case 3, 
        else:
            bacteria_dict[pathway] = {"Index": i, "Common-Name": row['Common-Name'], "Common-Name2": row['Common-Name2'], "Classification": "Unclassified"}
    
    # Sanity check
    i = 0
    j = 0
    k = 0
    
    for pathway,data in bacteria_dict.items():
        if data['Classification'] == "Unclassified":
             i = i+1
        elif data['Classification'] == "Bacteria":
            j = j+1
        elif data['Classification'] == "Non-Bacteria":
            k = k+1
        else:
            print("error in classification")
    
    print(i,j,k)


    df_bacterial = pd.DataFrame(bacteria_dict).T
    df_bacterial = df_bacterial.reset_index()
    df_bacterial = df_bacterial.rename(columns={"index":"Pathways"})
    df_bacterial = df_bacterial.set_index('Index')
    #df_bacterial.columns = ['Pathways_Index', 'Common-Name', 'Common-Name2', 'Classification'] 
    df_bacterial.reset_index()
    return(df_bacterial)
    

def merging_files_pathways(df_pwy_bacteria, pwy_file1, pwy_file2):
    df_file1 = pd.read_csv(pwy_file1, sep="\t")
    df_file2 = pd.read_csv(pwy_file2, sep="\t")
    df_columns_file1 = ["Reaction-List","EC-Number","Ontology - pathway type"]
    for column in df_columns_file1:
        df_file1[column] = df_file1[column].str.strip()
        df_file1[column] = df_file1[column].str.replace(' // ',';')
    
    df_merged = pd.merge(df_pwy_bacteria, df_file1, how = "inner", on = "Pathways")
    df_merged = pd.merge(df_merged, df_file2, how = "inner", on = "Pathways")
    
    # I only need to add the last file that has the summary
    
    
    
    
    
# %%
def main():
    df = filtering_by_bacteria()
    merging_files_pathways(df, "pathways_reactions_ECs_category.txt","Names_and_description_pwys.txt")
    #print(df.head(10))
    


# %%
main()
# %%
