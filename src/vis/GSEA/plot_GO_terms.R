#!/usr/bin/env Rscript
suppressPackageStartupMessages({
   library(here)
   library(tidyverse)
   library(topGO)
   library(Rgraphviz)
   })

snakemake@source(here("src","utils","logging.R"))
snakemake@source(here("src","utils","topGO.R"))
snakemake@source(here("src","utils","fig_params.R"))


geneID2GO <- readMappings(file = snakemake@input[["geneGOmappings"]])

sig_genes_df <- read_tsv(snakemake@input[["sigGenelist"]],
                        col_types = "cdddddd") %>% 
                        dplyr::select(Genes,padj)

print(head(sig_genes_df))

topGODataObject <- new("topGOdata", 
                       ontology = snakemake@wildcards[["GO_category"]], 
                       description = "DE genes",
                       allGenes = getGeneList(sig_genes_df,geneID2GO),
                       annot = annFUN.gene2GO,
                       gene2GO = geneID2GO,
                       nodeSize = as.integer(snakemake@params[["nodeSize"]]))

resultFisher <- runTest(object = topGODataObject,
                        algorithm = snakemake@wildcards[["algo"]], 
                        statistic = "fisher")

svg(snakemake@output[["GO_graph"]], 
    width = width_in,
    height = height_in)

showSigOfNodes(GOdata = topGODataObject, 
               score(resultFisher), 
               firstSigNodes = as.integer(snakemake@params[["firstSigNodes"]]), 
               useInfo ='all')
dev.off()
