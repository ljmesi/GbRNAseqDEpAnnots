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

replace_row_names_with_col <- function(df, col_name_to_row_name){
  df %<>% remove_rownames %>% 
        column_to_rownames(var = col_name_to_row_name)
}

create_meta <- function(df){
  sample_names <- paste(df$Condition, 
                        df$Sample_ID, 
                        sep = "-")
  condition <- paste(df$Condition, 
                     df$Phenotype, 
                     sep = "_")
  metadata <- data.frame(sample_names, condition, 
                         row.names = 1, 
                         stringsAsFactors = TRUE)
}

# Form a sample names vector
get_sample_names <- function(metadf){
    sample_names <- paste(meta$Condition, 
                        meta$Sample_ID, 
                        sep = ": ")
}


