#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

snakemake@source(here("src","utils","fig_params.R"))

### Set a colour palette
heat_colours <- brewer.pal(9, "Blues")

dds <- readRDS(snakemake@params[[1]])

transformed_data <- readRDS(snakemake@input[[1]])

# Indexes of 20 highest means of normalised counts of all samples from all genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

# Tell which samples belong to which condition (important for annotation)
condition <- colData(dds)[,c("condition")]
annotation <- as.data.frame(condition)
rownames(annotation) <- rownames(colData(dds))

# Open file
svg(snakemake@output[["count_matrix"]], 
    width = width_in, 
    height = height_in)


normalised_counts <- assay(transformed_data)


# Create the plot
pheatmap(normalised_counts[select,],
         color = heat_colours,
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=annotation)

# Close the file
dev.off()

# Calculate euclidian distance matrix between all samples
sampleDists <- dist(t(normalised_counts), method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- row.names(sampleDistMatrix)
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

# Open file
svg(snakemake@output[["sample_to_sample_distances"]], 
    width = width_in, 
    height = height_in)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Close the file
dev.off()

print("DESeq2 model dispersion plot drawn")
