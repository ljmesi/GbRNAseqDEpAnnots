#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

snakemake@source(here("src","utils","fig_params.R"))


vsd <- readRDS(snakemake@input[[1]])

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#### Plot sample to sample distance matrix
svg(snakemake@output[[1]], 
    width = width_in, 
    height = height_in)


pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         show_colnames = TRUE,
         show_rownames = TRUE)

dev.off()
