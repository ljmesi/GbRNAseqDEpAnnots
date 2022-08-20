.ONESHELL:
SHELL = /bin/bash
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: \
run \
clean \
report \
summary \
help \
update_env

# The used conda env definition file is located in "envs/snakemake.yml"
CURRENT_CONDA_ENV_NAME = snakemake
ACTIVATE_CONDA = source $$(conda info --base)/etc/profile.d/conda.sh
CONDA_ACTIVATE = $(ACTIVATE_CONDA) ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)


ARGS = --forceall

## run: Run the snakemake pipeline
run:
	$(CONDA_ACTIVATE)
	rm -f data/interm/annotations/geneIDs_GOs.tsv
	snakemake \
	--use-conda \
	--cores 7 \
	$(ARGS)

## clean: Remove output files
clean:
	$(CONDA_ACTIVATE)
	snakemake --cores 1 --delete-all-output --verbose
	cp logs/README.md log_readme.md
	cp bmarks/README.md bmark.readme.md
	rm --recursive --force --verbose \
	logs/* \
	bmarks/*
	mkdir -p logs bmarks
	mv log_readme.md logs/README.md
	mv bmark.readme.md bmarks/README.md


## report: Make automatic snakemake report
report:
	$(CONDA_ACTIVATE)
	rm -f report.html
	snakemake --cores 1 \
	--report report.html

## summary: Make a summary
summary:
	$(CONDA_ACTIVATE)
	snakemake --summary

## help: Show this message
help:
	@grep '^##' ./Makefile

## update_env: Update conda environment to the latest version defined by env.yml file
update_env:
	$(ACTIVATE_CONDA)
	mamba env update --file env.yml
