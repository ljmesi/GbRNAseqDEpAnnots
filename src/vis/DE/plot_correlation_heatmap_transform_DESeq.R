#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(vsn)

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
                             annotation = select(mdata, condition))
    dev.off()
}

plot_meanSd <- function(out_file, 
                        transformed_obj, 
                        w_in = width_in, 
                        h_in = height_in){
    # Open file
    svg(out_file, 
        width = w_in, 
        height = h_in)

    # Create the plot
    meanSdPlot(assay(transformed_obj))
    # Close the file
    dev.off()
}

#### Read in meta data ####
meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")

metadata <- create_meta(meta)

dds <- readRDS(snakemake@input[["normalised_counts_DESeq"]])


#### Transform the data so that it's better visualised with heatmaps ####
rld <- rlog(dds, blind=TRUE)
vstd <- vst(dds, blind=TRUE)
# this gives log2(n + 1)
ntd <- normTransform(dds)

#### Store the transformed objects for further use ####
saveRDS(rld, file = snakemake@output[["rlog_transformed"]])
saveRDS(vstd, file = snakemake@output[["vs_transformed"]])
saveRDS(ntd, file = snakemake@output[["norm_transformed"]])

#### Plot meanSd ####
# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#effects-of-transformations-on-the-variance
plot_meanSd(snakemake@output[["rlog_meanSdPlot"]], rld)
plot_meanSd(snakemake@output[["vs_meanSdPlot"]], vstd)
plot_meanSd(snakemake@output[["norm_meanSdPlot"]], ntd)

#### Do the actual plotting ####
plot_heatmap(mat_cor = compute_pairwise_correlation(rld),
             output_dest = snakemake@output[["rlog_transformed_corr_heatmap"]],
             mdata = metadata)
plot_heatmap(mat_cor = compute_pairwise_correlation(rld),
             output_dest = snakemake@output[["vs_transformed_corr_heatmap"]],
             mdata = metadata)
