#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)

snakemake@source(here("src","utils","utils.R"))

padj_limit <- snakemake@params[["p_adj_limit"]]
print("Adjusted p-value limit used:")
padj_limit

# Prepare a dataframe with outliers (indicated by having a pvalue)
dds <- readRDS(snakemake@input[["DESeq_analysis_obj"]])

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

contrast <- c("Condition", snakemake@params[["contrast"]])

res_w_outliers <- results(dds, 
                          contrast = contrast, 
                          alpha = padj_limit,
                          cooksCutoff = FALSE, 
                          parallel = parallel,
                          pAdjustMethod = snakemake@params[["pAdjustMethod"]])

rowname_col <- "Genes"

res_w_outliers <- as.data.frame(res_w_outliers[order(res_w_outliers$padj),]) %>%
        rownames_to_column(rowname_col)

# Prepare a dataframe without outliers (indicated by having pvalue equal to NA)

res_without_outliers <- as.data.frame(readRDS(snakemake@input[["DESeq_results_shrinked"]])) %>%
        rownames_to_column(rowname_col)



#### Read in raw counts to a data frame ####
raw_counts <- read.table(file = snakemake@input[["raw_counts"]],
                header = TRUE,
                sep = "\t",
                row.names = 1)

meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")


colnames(raw_counts) <- c(rowname_col, get_sample_names(meta))

# Find rows which differ on pvalue column
outliers <- anti_join(x = res_w_outliers, 
                      y = res_without_outliers,
                      by = "pvalue"
                      ) %>%
# Filter to only rows which exist in both dataframes
        inner_join(raw_counts, by = rowname_col)

write.table(x = outliers, 
            file = snakemake@output[["table_outliers"]],
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
