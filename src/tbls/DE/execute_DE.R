#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)

dds <- readRDS(snakemake@input[["DESeq_raw"]])

padj_limit <- snakemake@params[["p_adj_limit"]]
print("Adjusted p-value limit used:")
padj_limit

# From https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/scripts/deseq2-init.R
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# dds <- dds[ rowSums(counts(dds)) >= 10, ]
# remove uninformative rows with maximum of 1 counts across all samples
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds_run <- DESeq(dds, parallel=parallel)
# dds_run <- DESeq(dds, sfType = "poscounts", parallel=parallel)

saveRDS(dds_run, file = snakemake@output[["DESeq_analysis_obj"]])

print("DESeq analysis run and the output saved to .RDS object") 


#### Non-shrinked log2 fold changes ####

# Obtain results as log2fold changes relative to the 12:12 condition
print("Print contrast parameters:")
print(snakemake@params[["contrast"]])

contrast <- c("condition", snakemake@params[["contrast"]])
print(paste0("This is the numerator: ", contrast[2], " and this is the denominator: ", contrast[3]))

DESeq_res <- results(dds_run, 
                     contrast = contrast, 
                     alpha = padj_limit, 
                     parallel = parallel,
                     pAdjustMethod = snakemake@params[["pAdjustMethod"]])
# Order results by adjusted p-values
res <- DESeq_res[order(DESeq_res$padj),]

saveRDS(res, 
        file = snakemake@output[["DESeq_results_non_shrinked"]])

print("DESeq analysis result obtained and saved to .RDS format") 

rowname_col <- "Genes"

write.table(x = as.data.frame(res) %>% 
                rownames_to_column(var = rowname_col), 
            file = snakemake@output[["tbl_non_shrinked"]],
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")


# Take the subset of padj < padj_limit
res_filtered <- subset(res, padj < padj_limit)

write.table(x = as.data.frame(res_filtered) %>% 
                rownames_to_column(var = rowname_col), 
            file = snakemake@output[["tbl_non_shrinked_padj_filtered"]],
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")




#### Shrinked log2 fold changes ####
# https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/workflow/scripts/deseq2.R
# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial
# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
DESeq_res_shrinked <- lfcShrink(dds = dds_run,
                                contrast = contrast,
                                res = DESeq_res,
                                type = "normal",
                                parallel = parallel)

# sort by p-value
res <- DESeq_res_shrinked[order(DESeq_res_shrinked$padj),]

saveRDS(res, file = snakemake@output[["DESeq_results_shrinked"]])

write.table(x = as.data.frame(res) %>% 
                rownames_to_column(var = rowname_col), 
            file = snakemake@output[["tbl_shrinked"]],
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Take the subset of padj < padj_limit

res_filtered <- subset(res, padj < padj_limit)

write.table(x = as.data.frame(res_filtered) %>% 
                rownames_to_column(var = rowname_col), 
            file = snakemake@output[["tbl_shrinked_padj_filtered"]],
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

print("The shrinked log2fold DESeq results obtained and saved to .RDS object")

print("mcols-command: column metadata")
mcols(res, use.names = TRUE)
print("Standard output of the DESeqResults object")
res
print("summary-command: basic data about the results")
summary(res)