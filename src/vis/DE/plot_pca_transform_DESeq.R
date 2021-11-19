#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(ggrepel)

set.seed(1)

save_fig <- function(file_name, plot_object){
  ggsave(
    filename = file_name,
    plot = plot_object,
    device = "svg",
    width = width_in, 
    height = height_in,
    units = "in",
    dpi = dpi_no,
    limitsize = TRUE)
}

create_pca_plot <- function(transformed_obj, plot_title, samples){

  sample_labels <- data.frame(samples)

  # Plot PCA
  pca_plot <- plotPCA(transformed_obj, intgroup = "condition") +
                      labs(title = plot_title)

  pca_plot$data <- cbind(pca_plot$data,sample_labels)

  # # Add text lables to the points
  pca_plot <- pca_plot +
    geom_text_repel(aes(label = samples))
}

snakemake@source(here("src","utils","fig_params.R"))

meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")
samples = meta$Sample_ID
print(samples)


#### Use transformed data which works better with PCA ####
rld <- readRDS(snakemake@input[["rlog_transformed"]])
rld_plt <- create_pca_plot(rld, "PCA of normalised and sample-wise rlog transformed counts", samples)
save_fig(snakemake@output[["rlog_transformed_pca"]], rld_plt)


vstd <- readRDS(snakemake@input[["vs_transformed"]])
vstd_plt <- create_pca_plot(vstd, "PCA of normalised and sample-wise variance stabilisation transformed counts", samples)
save_fig(snakemake@output[["vs_transformed_pca"]], vstd_plt)
