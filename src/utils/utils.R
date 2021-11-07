#!/usr/bin/env Rscript

merge_files <- function(dir_path, regex, header, cols){
    file_list <- list.files(dir_path, recursive = FALSE, pattern = regex)
    for (file in file_list){
        file_path <- paste0(dir_path,
                            file)
        # if the merged dataset doesn't exist, create it
        if (!exists("dataset")){
            dataset <- read_tsv(file_path, 
                                col_names = header,
                                col_types = cols)
        }

        # if the merged dataset does exist, append to it
        if (exists("dataset")){
            temp_dataset <- read_tsv(file_path, 
                                    col_names = header,
                                    col_types = cols)
            print(paste0("Tsv file: ", file_path, " read. Glimpse of it:"))
            glimpse(temp_dataset)

            dataset<-bind_rows(dataset, temp_dataset)
            rm(temp_dataset)
        }
    }
    dataset
}

create_meta <- function(meta, sample_names){
    Condition <- paste(meta$Condition, 
                        meta$Phenotype, 
                        sep = " ")

    metadata <- data.frame(sample_names,Condition)

    metadata %<>% remove_rownames %>% 
        column_to_rownames(var = "sample_names")
}

# Form a sample names vector
get_sample_names <- function(metadf){
    sample_names <- paste(meta$Condition, 
                        meta$Sample_ID, 
                        sep = ": ")
}


