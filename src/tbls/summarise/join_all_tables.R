#!/usr/bin/env Rscript
suppressPackageStartupMessages({
   library(here)
   library(tidyverse)
   })

snakemake@source(here("src","utils","logging.R"))

join_col_name <- "Genes"

fileNames = c(snakemake@input[["BP"]], 
              snakemake@input[["CC"]], 
              snakemake@input[["MF"]])

tables <- lapply(X = fileNames, 
                 FUN = read_tsv, 
                 col_names = TRUE,
                 col_types = "ccciididddddd")

sig_genes <- read_tsv(snakemake@input[["SG"]], 
                      col_names = TRUE)

gff_annotations <- read_tsv(snakemake@input[["GFF"]], 
                            col_names = TRUE)

allData <- full_join(tables[[1]], tables[[2]], by = join_col_name) %>%
    full_join(tables[[3]], by = join_col_name) 

allData <- sig_genes %>%
    left_join(allData, by = join_col_name) %>%
    left_join(gff_annotations, by = join_col_name)

write.table(allData, 
            file = snakemake@output[[1]], 
            sep = "\t",
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
