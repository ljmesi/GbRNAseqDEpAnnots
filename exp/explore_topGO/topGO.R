testNumSig <- function(resultWeightObject){
  # see how many results we get where weight01 gives a P-value <= 0.05:
  summary_results_def <- summary(attributes(resultWeightObject)$score <= snakemake@config[["GSEA"]][["p_value"]])
  defaultNumSig <- 20
  # Initilise the variable
  numsignif <- defaultNumSig
    
  # If there were no significant GO-terms according to the default algorithm
  if (length(summary_results_def) == 2){
      numsignif <- defaultNumSig
    }else{
      # how many terms is it true that P <= 0.05
      numsignif <- as.integer(summary_results_def[[3]])
      # We'd like to have at least the defaultNumSig amount reported
      if (numsignif < defaultNumSig){
        numsignif <- defaultNumSig
      }
    }
  numsignif
}

getGeneList <- function(fileName,
                        GOmappings){
  # Define the genes of interest for topGO
  genesOfInterest <- read.table(fileName,
                                header=TRUE)
  genesOfInterest <- as.character(genesOfInterest$Genes)
  geneUniverse <- names(GOmappings)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  length(genesOfInterest)
  geneList
}


runTestWrap <- function(algo,topGODataO,stat) {
   return(runTest(object = topGODataO,
                  algorithm = algo, 
                  statistic = stat))
  }

runFisher <- function(algos = list("classic","elim","weight01","parentchild","weight","lea"), 
                      topGODataObject, 
                      stat = "fisher") {
  results <- lapply(X = algos, 
                    FUN = runTestWrap,
                    topGODataO = topGODataObject,
                    stat = stat)
}

genDF <- function(results, 
                  topGODataObject, 
                  goTermCharLen = 100){
  allRes <- GenTable(topGODataObject,
                     classicFisher = results[[1]],
                     elimFisher = results[[2]],
                     weight01Fisher = results[[3]],
                     parentchildFisher = results[[4]],
                     weightFisher = results[[5]],
                     leaFisher = results[[6]],
                     orderBy = "weight01Fisher",
                     ranksOf = "classicFisher",
                     numChar = goTermCharLen,
                     topNodes = testNumSig(results[[3]]))
}

