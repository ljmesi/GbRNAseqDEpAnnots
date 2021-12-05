#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)

snakemake@source(here("src","utils","fig_params.R"))

set.seed(1)

create_gene_symbols_list <- function(gene_annotations_df, names = "Genes", values = "Gene_name"){
  # Building a name/value association list
  IDs_symbols_mapping <- list()
  keys <- as.vector(gene_annotations_df[[names]])
  vals <- as.vector(gene_annotations_df[[values]])
  IDs_symbols_mapping[keys] <- vals
  return(IDs_symbols_mapping)
}

fill_nas_with_gene_IDs <- function(IDs_symbols_mapping){
  for (i in  1:length(IDs_symbols_mapping)){
    if(is.na(IDs_symbols_mapping[[i]])){
      gene_ID <- names(IDs_symbols_mapping)[i]
      IDs_symbols_mapping[i] <- gene_ID
    }
  }
  return(IDs_symbols_mapping)
}

find_matching_indexes <- function(all_geneIDs, de_geneIDs_mapping){
  gene_names <- unlist(de_geneIDs_mapping)
  gene_ID_indexes <- match(all_geneIDs, names(gene_names))
  return(gene_ID_indexes)
}


replace_gene_IDs_with_gene_symbols <- function(indexes_of_DE_genes, all_geneIDs, IDs_symbols_mapping){
  # Go through all indexes of DE genes in the main geneIDs vector
  for (i in  1:length(indexes_of_DE_genes)){
    index <- indexes_of_DE_genes[i]
    # Assign gene symbol (if it exists) in the place of gene ID
    if(!is.na(index)){
      all_geneIDs[i] <- IDs_symbols_mapping[[index]]
    }
  }
  return(all_geneIDs)
}


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

# Redo volcano plot
res <- readRDS(snakemake@input[["DESeq_results_shrinked"]])
row_names <- rownames(res)

gene_annotations_df <- read.table(
  file = snakemake@input[["gene_annotations"]], 
  sep = '\t', 
  header = TRUE, 
  stringsAsFactors = FALSE)

de_geneIDs_symbols_mapping <- fill_nas_with_gene_IDs(
  create_gene_symbols_list(
    gene_annotations_df))
gene_ID_indexes <- find_matching_indexes(row_names, de_geneIDs_symbols_mapping)
rownames(res) <- replace_gene_IDs_with_gene_symbols(gene_ID_indexes, row_names, de_geneIDs_symbols_mapping)

#### Plot sample to sample distance matrix
svg(snakemake@output[["volcano_plot_pub"]], 
    width = width_in, 
    height = height_in)

EnhancedVolcano(res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  ylim = c(0,12.7),
  xlim = c(-6,6),
  ylab = bquote(-~Log[10]~ 'Adj P'),
  title = "Volcanoplot",
  subtitle = "",
  legendPosition = 'right',
  legendLabels=c(' Not sig.',
                 bquote(~Log[2]~'FC'), 
                 ' p-value',
                 bquote(' p-value &' ~Log[2]~'FC')),
  legendLabSize = 12,
  col = c("#000000", "#E69F00", "#56B4E9", "#882255"),
  drawConnectors = TRUE,
  widthConnectors = 0.75
  )

dev.off()
