.ONESHELL:
SHELL = /bin/bash
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: run \
clean \
report \
help

# The conda env definition file "env.yml" is located in the project's root directory
CONDA_ENV_NAME = each_gene_into_one_row
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CONDA_ENV_NAME)

#ARGS = --forceall --rerun-incomplete
ARGS = --forceall

start_spyder:
	$(CONDA_ACTIVATE)
	spyder each_gene_in_one_row.py

script_run:
	$(CONDA_ACTIVATE)
	python \
	each_gene_in_one_row.py \
	combined-genes-GO_Fisher.tsv \
	test.tsv



## help: Show this message
help:
	@grep '^##' ./Makefile
