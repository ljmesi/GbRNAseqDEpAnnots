#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

snakemake@source(here("src","utils","logging.R"))

library(tidyverse)

snakemake@source(here("src","utils","fig_params.R"))

#### Read in raw counts to a data frame ####
raw_counts <- read.table(file = snakemake@input[["raw_counts"]],
                header = TRUE,
                sep = "\t",
                row.names = 1)

# Calculate mean for each gene (each row)
mean_counts <- apply(raw_counts[, 2:5], 1, mean)
# Calculate variance for each gene (each row)
variance_counts <- apply(raw_counts[, 2:5], 1, var)

# Create data frame with mean and variance for every gene
df <- data.frame(mean_counts, variance_counts)

plot <- ggplot(df) +
  geom_point(aes(x = mean_counts, 
                 y = variance_counts
                )
            ) +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "Log mean vs log variance of raw counts of condition 12/12",
       x = "Mean counts per gene",
       y = "Variance per gene"
      )

ggsave(
  filename = snakemake@output[["pl"]],
  plot = plot,
  device = "svg",
  width = width_in,
  height = height_in,
  units = "in", #c("in", "cm", "mm"),
  # dpi = dpi_no,
  limitsize = TRUE
)
