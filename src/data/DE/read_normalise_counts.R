#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(DESeq2)
library(tidyverse)
library(magrittr)

snakemake@source(here("src","utils","utils.R"))

#### Read in the raw data and experimental metadata ####
file_name <- snakemake@input[["raw_counts_table"]]
                 
meta <- read_tsv(file = snakemake@input[["experiment_metadata"]],
                 col_types = "ccc")

if(file_name == "data/raw/counts.raw.txt") {
    raw <- read_delim(file = file_name, 
                    delim = " ",
    # The second column is char here because of the issue with the counts table
                    col_types = "cciiiiiii")
    # Rename all columns to nicer ones
    colnames(raw) <- c("Gene", get_sample_names(meta))
    # Assign where erroneous " " existed in a column
    erroneous_col_index <- 2
    # Convert the column with issues (an extra space in front of the integer) to integer
    # Create a new column where the erroneous space existed
    raw %<>% mutate(despaced = gsub(" ", "", .[[erroneous_col_index]]))
    # Convert the despaced column to integer and assign to its old location
    raw$despaced %<>% as.integer
    raw[,erroneous_col_index] <- raw$despaced
    # Remove the helper column
    raw %<>% select(-despaced)
} else if(file_name == "data/raw/STAR_counts_all.csv") {
    raw <- read_delim(file = file_name, 
                    delim = ",",
                    col_types = "ciiiiiiii")
    # Rename all columns to nicer ones
    colnames(raw) <- c("Gene", get_sample_names(meta))
} else if(file_name == "data/raw/subread_counts.txt") {
    raw <- read_delim(file = file_name, 
                      delim = "\t",
                      comment = "#",
                      col_types = "ccccciiiiiiiii") %>%
        select(-c("Chr","End","Start","Strand","Length"))
    # Rename all columns to nicer ones
    colnames(raw) <- c("Gene", get_sample_names(meta))    
} else {
    stop("This input file is unknown")
}


# Log a glimpse of the counts data
glimpse(raw)


#### Store the raw data ####
# Save the raw data for use with finding outliers and 
# for plotting log mean log variance plot
write.table(x = raw, 
            file = snakemake@output[["raw_counts"]], 
            append = FALSE, 
            sep = "\t", 
            dec = ".",
            row.names = TRUE, 
            col.names = TRUE,
            quote = FALSE)


raw %<>% remove_rownames %>% 
    column_to_rownames(var = "Gene")

#### Create DESeq object ####
dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = create_meta(meta, 
                                                    paste(meta$Condition, 
                                                          meta$Sample_ID, 
                                                          sep = ": ")),
                              design = ~ Condition)

# Save the DESeq object for future use
saveRDS(dds, file = snakemake@output[["DESeq_obj"]])
print("Raw counts read to a DESeq2 object")


#### Check that raw counts and experiment design dataframes are correct ####
print(
    paste0("All row and columns exist in raw data and design data frame: ",
            all(colnames(raw) %in% rownames(metadata))
          )
     )
print(
    paste0("All row names are in same order as column names: ",
            all(colnames(raw) == rownames(metadata))
          )
     ) 


#### Prefilter rows with no information and normalise the data ####
# remove uninformative rows with maximum of 1 counts across all samples 
dds <- dds[rowSums(counts(dds)) > 1, ]

dds_normalised <- estimateSizeFactors(dds)

print("Counts normalised")


#### Store intermediary files ####
saveRDS(dds_normalised, 
        file = snakemake@output[["normalised_counts_DESeq"]])

write.table(x = counts(dds_normalised, 
                       normalized=TRUE), 
            file = snakemake@output[["normalised_counts"]], 
            append = FALSE, 
            sep = "\t", 
            dec = ".",
            row.names = TRUE, 
            col.names = TRUE,
            quote = FALSE)


#### Do some logging ####
print("Normalised counts stored in a tsv file")

print(paste0("Normalisation factors applied to each sample are: ",
             sizeFactors(dds_normalised)
            ))