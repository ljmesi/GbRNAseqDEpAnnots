.ONESHELL:
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: start, \
help

CURRENT_CONDA_ENV_NAME = download_GO_info
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

OUTFILE = all_samples_genes_count.tsv

## download_GO_info: Download GO info
download_GO_info:
	$(CONDA_ACTIVATE)
	find STAR_counts/ -name "*.tab" -print0 | xargs -0 ./collate_counts.py -o $(OUTFILE)

## help: Show this message
help:
	@grep '^##' ./Makefile
