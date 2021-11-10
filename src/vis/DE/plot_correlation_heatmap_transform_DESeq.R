#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

snakemake@source(here("src","utils","fig_params.R"))
snakemake@source(here("src","utils","utils.R"))

#### Set a colour palette ####
heat_colours <- brewer.pal(6, "YlOrRd")


#### Calculate pairwise correlations ####
compute_pairwise_correlation <- function(DESeq_object) {
    # Extract the transformed matrix from the object
    mat <- assay(DESeq_object)
    # Compute pairwise correlation values
    mat_cor <- cor(mat)
    return(mat_cor)
}


#### Plot heatmap ####
plot_heatmap <- function(mat_cor, 
                         output_dest,
                         colours = heat_colours, 
                         mdata = metadata, 
                         w_in = width_in,
                         h_in = height_in) {
    svg(output_dest, 
        width = w_in,
        height = h_in)

    # Plot heatmap
    corr_heatmap <- pheatmap(mat_cor,
                            color = colours,  
                            annotation = select(mdata, 
                                                Condition))
    dev.off()
}


#### Read in meta data ####
meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")

metadata <- create_meta(meta, get_sample_names(meta))

dds <- readRDS(snakemake@input[["normalised_counts_DESeq"]])


#### Transform the data so that it's better visualised with heatmaps ####
rld <- rlog(dds, 
            blind=TRUE)


vstd <- vst(dds, 
           blind=TRUE)


#### Do the actual plotting ####
plot_heatmap(mat_cor = compute_pairwise_correlation(rld),
             output_dest = snakemake@output[["rlog_transformed_corr_heatmap"]]
            )

plot_heatmap(mat_cor = compute_pairwise_correlation(rld),
             output_dest = snakemake@output[["vs_transformed_corr_heatmap"]]
            )
