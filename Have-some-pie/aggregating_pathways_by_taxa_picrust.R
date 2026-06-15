# Daniel Castaneda Mogollon, PhD
# May 20th, 2026
# This script aggregates the pathways for each taxa, instead of having them as ASV contribution.
# then, the output will serve for Maaslin2 to identify significantly enriched pathways at the taxa level and not community-based level.

library(dplyr)

path = "/Users/danielcm/Desktop/SickKids/PICRUSt2.6/"
communities = c("ns1","ns6","s2","s5")
setwd(path)

for (community in communities){
    in_file = read.csv(paste0(community, "_output/", community, "_pathway_inference/path_only_abun_strat_no_taxa.tsv"), sep = "\t", check.names = FALSE)
    #in_file$pathway_taxa = paste(in_file$pathway, in_file$Taxa, sep = "_")

    df = as.data.frame(in_file)

    #print(colnames(df))
    df_aggregated = df %>%
        #group_by(pathway, Taxa) %>%
        group_by(pathway) %>%
        summarize(
            across(where(is.numeric), \(x) sum(x, na.rm = TRUE))
        )

    # Sanity check, the sum of the original df should match the sum of the aggregated df
    #col_sums = colSums(df[,5:ncol(df)-1])
    colnames(df[,2:ncol(df)])
    col_sums = colSums(df[,2:ncol(df)])
    sum = sum(col_sums)


    #col_sums2 = colSums(df_aggregated[,3:ncol(df_aggregated)])
    #colnames(df_aggregated[,3:ncol(df_aggregated)])
    col_sums2 = colSums(df_aggregated[,2:ncol(df_aggregated)])
    colnames(df_aggregated[,2:ncol(df_aggregated)])

    sum2 = sum(col_sums2)
    sum2

    if(sum != sum2){
        print(paste0("The sums do not match, ",sum," vs ",sum2," there might be an error in the aggregation process."))
        print(sum2-sum)
    } else {
        print(paste0("The sums match, ",sum," vs ",sum2," the aggregation process is likely correct."))
        print(sum2-sum)
    }
    # Sanity check passed

    write.table(df_aggregated, paste0(community, "_output/", community, "_pathway_inference/path_only_abun_strat_no_taxa_aggregated.tsv"), sep = "\t", row.names = FALSE)
}
