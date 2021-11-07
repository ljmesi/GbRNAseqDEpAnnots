#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)

snakemake@source(here("src","utils","fig_params.R"))

dds_run <- readRDS(snakemake@input[["DESeq_results"]])

# Open file
svg(snakemake@output[["MAplot"]], 
    width = width_in, 
    height = height_in
   )
# Create the plot
plotMA(dds_run)
# Close the file
dev.off()

print(
    paste0("MAplot with ",
           snakemake@wildcards[["shrink_status"]],
           " DE log2fold values drawn"
          )
    )
