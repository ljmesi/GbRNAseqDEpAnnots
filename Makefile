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
CONDA_ENV_NAME = snakemake
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CONDA_ENV_NAME)

#ARGS = --rerun-incomplete

## run: Run the snakemake pipeline
run:
	$(CONDA_ACTIVATE)
	snakemake $(ARGS) \
	--cores 7 \
	--use-conda 

## clean: Remove output files
clean:
	$(CONDA_ACTIVATE)
	snakemake --cores 1 --delete-all-output
	rm -f \
	"data/interm/DE/normalised_counts.tsv" \
	"results/figs/DE/normalised.correlation.rlog.svg" \
	"results/tbls/DE/shrinked_padj-filtered.tsv" \
	"results/figs/DE/log10_cooks_distances.svg" \
	"results/figs/DE/deseq_dispersions.svg" \
	"results/tbls/DE/outliers.tsv", \
	"results/figs/DE/MAplot_shrinked.svg", \
	"results/figs/DE/volcanoplot.svg", \
	"results/figs/DE/DE_heatmap.svg", \
	"data/interm/DE/genelist_padj-filtered.tsv", \
	"data/interm/GSEA/geneIDs_GOs.tsv"

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
