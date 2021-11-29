#!/usr/bin/env Rscript
suppressPackageStartupMessages({
   library(here)
   library(tidyverse)
   library(topGO)
   library(magrittr)
   })

snakemake@source(here("src","utils","logging.R"))
snakemake@source(here("src","utils","topGO.R"))


#### Write GO terms ####

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

results <- runFisher(topGODataObject = topGODataObject)

resDF <- genDF(results = results,
               topGODataObject = topGODataObject)


#### Assign GO terms to GeneIDs ####

goTermsOfInterest <- pull(resDF,GO.ID)
numSignGOterms <- length(goTermsOfInterest)

go2geneIDs <- genesInTerm(topGODataObject, goTermsOfInterest)

GOtermGenes <- data.frame(GO.ID = character(numSignGOterms),
                          Genes = character(numSignGOterms),
                          stringsAsFactors = FALSE)

# From a really helpful blogpost:
# https://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
for (i in 1:numSignGOterms){
   goTerm <- goTermsOfInterest[i]
   GOtermGenes[i,1] <- goTerm
   mygenesforterm <- go2geneIDs[goTerm][[1]]
   # Change char vector to one continuous char
   mygenesforterm <- paste(mygenesforterm, collapse=',')
   GOtermGenes[i,2] <- mygenesforterm
}
# Merge GO term - Gene ID data.frame with the data with p values 
# and GO term descriptions, etc.
allData <- merge(resDF, GOtermGenes, by = "GO.ID")

write.table(allData, 
            file = snakemake@output[[1]], 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
