#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

snakemake@source(here("src","utils","utils.R"))
snakemake@source(here("src","utils","fig_params.R"))

padj_limit <- snakemake@params[["p_adj_limit"]]
print("Adjusted p-value limit used:")
padj_limit

### Set a colour palette
heat_colours <- brewer.pal(6, "YlOrRd")


#### Read in metadata ####
meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")

Condition <- paste(meta$Condition, 
                  meta$Phenotype, 
                  sep = " ")

# Form a sample name vector
sample_names <- paste(meta$Condition, 
                      meta$Sample_ID, 
                      sep = ": ")

annotation <- data.frame(sample_names,Condition)

annotation %<>% remove_rownames %>% column_to_rownames(var = "sample_names")


# Obtain the normalised count values of DE genes
results <- readRDS(snakemake@input[["DESeq_results_shrinked"]])

# The genelabels column should be now filled fully
# because it is needed in subsetting normalised_counts
results$genelabels <- rownames(results)

results <- as_tibble(results) %>% 
  column_to_rownames(var = "genelabels")

## Extract significant genes
sigOE <- subset(results, padj < padj_limit)

normalised_counts <- read.table(file = snakemake@input[["normalised_counts"]],
                                header = TRUE,
                                sep = "\t",
                                row.names = 1)


### Extract normalized expression for significant genes
norm_OEsig <- normalised_counts[rownames(sigOE),]

meta <- read_tsv(file = snakemake@input[["metadata"]],
                 col_types = "ccc")

colnames(norm_OEsig) <- get_sample_names(meta)

norm_OEsig <- as.matrix(norm_OEsig)

# Plot the DE heatmap
svg(snakemake@output[["DE_heatmap"]], 
    width = width_in, 
    height = height_in)

pheatmap(norm_OEsig,
         color = heat_colours, 
         cluster_rows = TRUE, 
         show_rownames =TRUE,
         annotation_col = annotation,
         fontsize_row = 5, 
         border_color=NA, 
         fontsize = 25,
         scale="row", 
         height = 20)
dev.off()
