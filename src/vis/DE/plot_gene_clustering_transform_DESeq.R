#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(vsn)
library(genefilter)

snakemake@source(here("src","utils","fig_params.R"))
snakemake@source(here("src","utils","utils.R"))

extract_n_highest_var_across_samples <- function(dESeqTransformObject, n){
  topVarGenes <- head(order(rowVars(assay(dESeqTransformObject)), decreasing = TRUE), n)
}

get_deviation_from_avg_across_samples <- function(dESeqTransformObject, pos_of_genes_of_interest){
  mat <-assay(dESeqTransformObject)[ pos_of_genes_of_interest, ]
  mat <- mat - rowMeans(mat)
}

#### Plot heatmap ####
plot_heatmap <- function(mat, 
                         output_dest,
                         annot = annotation,
                         heatmap_title = "Clustering of 20 genes with the variance across samples",
                         w_in = width_in,
                         h_in = height_in) {
  svg(output_dest, 
      width = w_in,
      height = h_in)
  
  # Plot heatmap
  corr_heatmap <- pheatmap(mat,
                         main = heatmap_title,
                         annotation = annot)
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

#### Read in data ####
dds <- readRDS(snakemake@input[["DESeq_counts"]])

#### Transform the data so that it's better visualised with heatmaps ####
rld <- rlog(dds, blind=TRUE)
vstd <- vst(dds, blind=TRUE)


#### Transform the data so that it's better visualised with heatmaps ####
rld <- rlog(dds, blind=TRUE)
vstd <- vst(dds, blind=TRUE)


#### Store the transformed objects for further use ####
saveRDS(rld, file = snakemake@output[["rlog_transformed"]])
saveRDS(vstd, file = snakemake@output[["vs_transformed"]])


#### Plot meanSd ####
# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#effects-of-transformations-on-the-variance
plot_meanSd(snakemake@output[["rlog_meanSdPlot"]], rld)
plot_meanSd(snakemake@output[["vs_meanSdPlot"]], vstd)


#### Do the actual plotting ####
#### Plot heatmap ####
plot_heatmap(mat = get_deviation_from_avg_across_samples(dESeqTransformObject = rld,
                                                         pos_of_genes_of_interest = extract_n_highest_var_across_samples(dESeqTransformObject = rld,
                                                                                                                         n = 20)),
             output_dest = snakemake@output[["rlog_gene_clustering_heatmap"]],
             annot = create_annotation(dESeqTransformObject = rld))

plot_heatmap(mat = get_deviation_from_avg_across_samples(dESeqTransformObject = vstd,
                                                         pos_of_genes_of_interest = extract_n_highest_var_across_samples(dESeqTransformObject = vstd,
                                                                                                                         n = 20)),
             output_dest = snakemake@output[["vs_gene_clustering_heatmap"]],
             annot = create_annotation(dESeqTransformObject = vstd))
