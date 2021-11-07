SHELL = /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# The conda env definition file "env.yml" is located in the project's root directory
CONDA_ENV_NAME = snakemake
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CONDA_ENV_NAME)

.PHONY: run, \
clean

run: clean
	$(CONDA_ACTIVATE) ; \
	snakemake --cores 7 --use-conda

clean:
	rm -f \
	"data/interm/DE/normalised_counts.tsv" \
	"results/figs/DE/normalised.correlation.rlog.svg" \
	"results/tbls/DE/shrinked_padj-filtered.tsv"
