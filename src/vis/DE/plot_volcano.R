#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(EnhancedVolcano)
})

snakemake@source(here("src","utils","fig_params.R"))

set.seed(1)



padj_limit <- snakemake@params[["p_adj_limit"]]
print(cat("Adjusted p-value limit used: ", padj_limit))
results <- readRDS(snakemake@input[["DESeq_results_shrinked"]])

gene_annotations_df <- read.table(
  file = snakemake@input[["gene_annotations"]], 
  sep = '\t', 
  header = TRUE, 
  stringsAsFactors = FALSE)



# Add a column labeling whether the gene is significant
#threshold_OE <- results$padj < padj_limit 
results$threshold <- results$padj < padj_limit
# Count how many TRUE values were found
topGeneNumber <- sum(results$threshold, na.rm = TRUE)
print(cat("How many are there which pass the adjusted p-value threshold of ", padj_limit, ":"))
print(topGeneNumber)

res_ordered <- results[order(results$padj), ]

## Create a column to indicate which genes to label
res_ordered$genelabels <- ""
res_ordered$genelabels[1:topGeneNumber] <- 
  rownames(res_ordered)[1:topGeneNumber]


# Plot volcano
res_ordered <- as_tibble(res_ordered)
# Remove NAs from padj column
res_ordered <- drop_na(res_ordered, padj)
plot <- ggplot(res_ordered) +
  geom_point(aes(x = log2FoldChange, 
                  y = -log10(padj), 
                  colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, 
                      y = -log10(padj), 
                      label = genelabels),
                  max.overlaps = Inf) +
  ggtitle("Overexpression of genes") +
  xlab("log2 fold change") + 
  ylab("-log10 transformed adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(plot_title_rel_size), 
                                  hjust = plot_title_justification),
        axis.title = element_text(size = rel(axis_title_rel_size))) 

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

print("Volcano plotted")

# Redo volcano plot
res <- readRDS(snakemake@input[["DESeq_results_shrinked"]])
row_names <- rownames(res)


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
