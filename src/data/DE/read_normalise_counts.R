#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(magrittr)

snakemake@source(here("src","utils","utils.R"))

#### Read in the raw data and experimental metadata ####
meta <- read_tsv(file = snakemake@input[["experiment_metadata"]],
                 col_types = "ccc")

raw <- read_delim(file = snakemake@input[["raw_counts_table"]], 
                delim = "\t",
                col_types = "ciiiiiiii")

print("Glimpse of the counts data:")
glimpse(raw)

#### Create DESeq object ####
dds <- DESeqDataSetFromMatrix(countData = replace_row_names_with_col(raw, "Genes"),
                              colData = create_meta(meta),
                              design = ~ condition)

print("List the coefficients:")
resultsNames(dds)

# Save the DESeq object for future use
saveRDS(dds, file = snakemake@output[["DESeq_obj"]])
print("Raw counts read to a DESeq2 object")


#### Check that raw counts and experiment design dataframes are correct ####
print(paste0("All row and columns exist in raw data and design data frame: ",
             all(colnames(raw) %in% rownames(metadata))))
print(paste0("All row names are in same order as column names: ",
             all(colnames(raw) == rownames(metadata)))) 


#### Prefilter rows with no information and normalise the data ####
# remove uninformative rows with maximum of 1 counts across all samples 
dds <- dds[rowSums(counts(dds)) > 1, ]

dds_normalised <- estimateSizeFactors(dds)

print("Counts normalised")


#### Store intermediary normalised counts for preview ####
write.table(x = counts(dds_normalised, 
                       normalized=TRUE), 
            file = snakemake@output[["normalised_counts"]], 
            append = FALSE, 
            sep = "\t", 
            dec = ".",
            row.names = TRUE, 
            col.names = TRUE,
            quote = FALSE)


#### Do some logging ####
print("Normalised counts stored in a tsv file")

print(paste0("Normalisation factors applied to each sample are: ",
             sizeFactors(dds_normalised)))
