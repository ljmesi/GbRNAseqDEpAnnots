#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(ggrepel)

snakemake@source(here("src","utils","fig_params.R"))

set.seed(1)

padj_limit <- snakemake@params[["p_adj_limit"]]
print("Adjusted p-value limit used:")
padj_limit

results <- readRDS(snakemake@input[["DESeq_results_shrinked"]])

topGeneNumber <- snakemake@params[["volcano_top_genes"]]

# Add a column labeling whether the gene is significant
#threshold_OE <- results$padj < padj_limit 
results$threshold <- results$padj < padj_limit 

res_ordered <- results[order(results$padj), ]

## Create a column to indicate which genes to label
res_ordered$genelabels <- ""
res_ordered$genelabels[1:topGeneNumber] <- 
  rownames(res_ordered)[1:topGeneNumber]


# Plot volcano
res_ordered <- as_tibble(res_ordered)
plot <- ggplot(res_ordered) +
  geom_point(aes(x = log2FoldChange, 
                 y = -log10(padj), 
                 colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, 
                      y = -log10(padj), 
                      label = genelabels)) +
  ggtitle("Overexpression of genes") +
  xlab("log2 fold change") + 
  ylab("-log10 transformed adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(plot_title_rel_size), 
                                  hjust = plot_title_justification),
        axis.title = element_text(size = rel(axis_title_rel_size))) 

print("Volcano plotted")

ggsave(
  filename = snakemake@output[["volcano_plot"]],
  plot = plot,
  device = "svg",
  width = width_in,
  height = height_in,
  units = "in", #c("in", "cm", "mm"),
  # dpi = dpi_no,
  limitsize = TRUE
)
