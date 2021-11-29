#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)

snakemake@source(here("src","utils","fig_params.R"))

dds_run <- readRDS(snakemake@input[["DESeq_analysis_obj"]])

# Open file
svg(snakemake@output[[1]], 
    width = width_in, 
    height = height_in)

# Create the plot
plotDispEsts(dds_run)
# Close the file
dev.off()

print("DESeq2 model dispersion plot drawn")
