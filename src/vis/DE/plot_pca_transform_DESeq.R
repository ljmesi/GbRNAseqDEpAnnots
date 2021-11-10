#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(ggrepel)

snakemake@source(here("src","utils","fig_params.R"))

set.seed(1)

dds <- readRDS(snakemake@input[["normalised_counts_DESeq"]])


meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")

sample_labels <- data.frame(samples = meta$Sample_ID)

print(sample_labels)

#### Transform the data so that it's better visualised with heatmaps ####
rld <- rlog(dds, 
            blind=TRUE)

# Plot PCA
pca_plot <- plotPCA(rld, intgroup = "Condition") +
                    labs(title = "PCA of normalised and sample-wise rlog transformed counts")

pca_plot$data <- cbind(pca_plot$data,sample_labels)

# # Add text lables to the points
pca_plot <- pca_plot +
  geom_text_repel(aes(label = samples))

ggsave(
  filename = snakemake@output[["rlog_transformed_pca"]],
  plot = last_plot(),
  device = "svg",
  width = width_in, 
  height = height_in,
  units = "in",
  dpi = dpi_no,
  limitsize = TRUE)


vstd <- vst(dds, 
           blind=TRUE)


pca_plot <- plotPCA(vstd, intgroup = "Condition") +
                    labs(title = "PCA of normalised and sample-wise variance stabilisation transformed counts")

pca_plot$data <- cbind(pca_plot$data,sample_labels)

# # Add text lables to the points
pca_plot <- pca_plot +
  geom_text_repel(aes(label = samples))

ggsave(
  filename = snakemake@output[["vs_transformed_pca"]],
  plot = last_plot(),
  device = "svg",
  width = width_in, 
  height = height_in,
  units = "in",
  dpi = dpi_no,
  limitsize = TRUE)
