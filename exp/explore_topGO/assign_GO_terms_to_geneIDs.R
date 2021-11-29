log <- file(snakemake@log[[1]], open="wt")
# Divert output to log file
sink(log)
# Divert messages, warning and stop to log file
sink(log, type="message")

library(tidyverse)
library(topGO)

snakemake@source("../../utils/topGO.R")

#### Write GO terms ####

geneID2GO <- readMappings(file = snakemake@input[["geneGOmappings"]])

topGODataObject <- new("topGOdata", 
                       ontology = snakemake@wildcards[["GO_category"]], 
                       description = "DE genes",
                       allGenes = getGeneList(snakemake@input[["sigGenelist"]],
                                              geneID2GO),
                       annot = annFUN.gene2GO,
                       gene2GO = geneID2GO,
                       nodeSize = 5)

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
            file = snakemake@output[["tsv_table"]], 
            sep = "\t",
            qmethod = "double",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)