# Daniel Castaneda Mogollon, PhD
# May 20th, 2026
# This script aggregates PICRUSt2 pathway abundances from ASV/sequence-level
# contributions to pathway-by-taxon contributions.

library(dplyr)

path <- "/Users/danielcm/Desktop/SickKids/PICRUSt2.6/"
communities <- c("ns1", "ns6", "s2", "s5")
setwd(path)

for (community in communities) {
    input_file <- paste0(community, "_output/", community, "_pathway_inference/path_abun_strat_with_taxa.tsv")
    output_file <- paste0(community, "_output/", community, "_pathway_inference/path_abun_strat_with_taxa_aggregated.tsv")

    in_file <- read.delim(input_file, check.names = FALSE)

    required_columns <- c("pathway", "sequence", "Taxa")
    missing_columns <- setdiff(required_columns, colnames(in_file))
    if (length(missing_columns) > 0) {
        stop(paste0(
            "Missing required column(s) in ", input_file, ": ",
            paste(missing_columns, collapse = ", ")
        ))
    }

    df_aggregated <- in_file %>%
        select(-sequence) %>%
        group_by(pathway, Taxa) %>%
        summarize(
            across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
            .groups = "drop"
        )

    sample_columns <- setdiff(colnames(in_file), required_columns)
    input_sum <- sum(colSums(in_file[, sample_columns, drop = FALSE]))
    aggregated_sum <- sum(colSums(df_aggregated[, sample_columns, drop = FALSE]))

    if (!isTRUE(all.equal(input_sum, aggregated_sum, tolerance = 1e-8))) {
        stop(paste0(
            "The sums do not match for ", community, ": ",
            input_sum, " vs ", aggregated_sum,
            ". There might be an error in the aggregation process."
        ))
    }

    print(paste0(
        community, ": aggregated ", nrow(in_file),
        " ASV-level rows into ", nrow(df_aggregated),
        " pathway-taxon rows; total abundance preserved."
    ))

    write.table(
        df_aggregated,
        output_file,
        sep = "	",
        row.names = FALSE,
        quote = FALSE
    )
}
