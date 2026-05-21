# Daniel Castaneda Mogollon, PhD
# May 20th, 2026
# This script aggregates the pathways for each taxa, instead of having them as ASV contribution.
# then, the output will serve for Maaslin2 to identify significantly enriched pathways at the taxa level and not community-based level.

library(dplyr)

path = "/Users/danielcm/Desktop/SickKids/PICRUSt2.6/"
communities = c("ns1","ns6","s2","s5")
setwd(path)

in_file = read.csv("ns1_output/ns1_pathway_inference/path_abun_strat_with_taxa.tsv", sep = "\t")
in_file$pathway_taxa = paste(in_file$pathway, in_file$Taxa, sep = "_")

df = as.data.frame(in_file)
df_aggregated = df %>%
    group_by(pathway, Taxa) %>%
    summarize(
        across(where(is.numeric), sum, na.rm = TRUE)
    )

# Sanity check, the sum of the original df should match the sum of the aggregated df
col_sums = colSums(df[,5:ncol(df)-1])
sum = sum(col_sums)
sum


col_sums2 = colSums(df_aggregated[,3:ncol(df_aggregated)])
colnames(df_aggregated[,3:ncol(df_aggregated)])

sum2 = sum(col_sums2)
sum2

if(sum != sum2){
    print("The sums do not match, there might be an error in the aggregation process.")
} else {
    print("The sums match, the aggregation process is likely correct.")
}
# Sanity check passed

write.table(df_aggregated, "ns1_output/ns1_pathway_inference/path_abun_strat_with_taxa_aggregated.tsv", sep = "\t", row.names = FALSE)
