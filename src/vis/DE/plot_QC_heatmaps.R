#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

snakemake@source(here("src","utils","fig_params.R"))
snakemake@source(here("src","utils","utils.R"))

### Set a colour palette
heat_colours <- brewer.pal(9, "Blues")

transformed_data <- readRDS(snakemake@input[[1]])
transformed_assay <- assay(transformed_data)

print("Transformed data assay:")
head(transformed_assay)

# Tell which samples belong to which condition (important for annotation)
annotation <- create_annotation(transformed_data)

# Indexes of "params" highest means of transformed counts of all samples from all genes
num_highest_means <- as.integer(snakemake@params[["num_highest_means"]])
select <- order(rowMeans(transformed_assay),decreasing=TRUE)[1:num_highest_means]

#### Plot count matrix ####
svg(snakemake@output[["count_matrix"]], 
    width = width_in, 
    height = height_in)

pheatmap(transformed_assay[select,],
         color = heat_colours,
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=annotation)
dev.off()

# Calculate euclidian distance matrix between all samples
sampleDists <- dist(t(transformed_assay), method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- row.names(sampleDistMatrix)
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

#### Plot sample to sample distance matrix
svg(snakemake@output[["sample_to_sample_distances"]], 
    width = width_in, 
    height = height_in)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

print("DESeq2 model dispersion plot drawn")
