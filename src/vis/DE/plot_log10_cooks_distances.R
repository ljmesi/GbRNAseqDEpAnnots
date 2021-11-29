#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)

snakemake@source(here("src","utils","fig_params.R"))

dds <- readRDS(snakemake@input[["DESeq_analysis_obj"]])

# Create data frame for previewing the cooks distances
df <- data.frame(assays(dds)[["cooks"]])
print(head(df))
glimpse(df)

# Open file
svg(snakemake@output[[1]], 
    width = width_in, 
    height = height_in)

# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#approach-to-count-outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

dev.off()