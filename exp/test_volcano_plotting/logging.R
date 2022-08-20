#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
# Divert output to log file
sink(log)
# Divert messages, warning and stop to log file
sink(log, type="message")
