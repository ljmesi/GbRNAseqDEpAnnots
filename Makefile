.ONESHELL:
SHELL = /bin/bash
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
	"results/tbls/DE/shrinked_padj-filtered.tsv" \
	"results/figs/DE/log_mean_log_variance.svg" \
	"results/figs/DE/deseq_dispersions.svg" \
	"results/tbls/DE/outliers.tsv", \
	"results/figs/DE/MAplot_shrinked.svg", \
	"results/figs/DE/volcanoplot.svg", \
	"results/figs/DE/DE_heatmap.svg", \
	"data/interm/DE/genelist_padj-filtered.tsv"
