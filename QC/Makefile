.ONESHELL:
SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules


#### NOTE: These rules in this script didn't work ####


QC_DIR = /crex/proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/QC

.PHONY: fastqc, \
seqkit , \
multiqc

fastqc:
	cd $(QC_DIR)
	module load bioinfo-tools gnuparallel/20180822 ; \
	cat ../data/sample_names.txt \
	| parallel --eta \
	--progress \
	--bar \
	--verbose \
	--jobs 8 \
	"singularity exec ../simgs/fastqc.sif fastqc ../data/raw/{}_R1_001.fastq.gz ../data/raw/{}_R2_001.fastq.gz --format fastq --outdir fastqc >> fastqc.log 2>&1"

seqkit:
	cd $(QC_DIR)
	module load bioinfo-tools gnuparallel/20180822 ; \
        cat ../data/sample_names.txt \
        | parallel --eta \
        --keep-order \
        --interactive \
        --progress \
        --bar \
        --verbose \
        --jobs 8 \
        -N1 \
        "singularity exec ../simgs/seqkit.sif seqkit stats ../data/raw/{}_R1_001.fastq.gz ../data/raw/{}_R2_001.fastq.gz >> seqkit_stats.tsv 2>&1"

multiqc:
	singularity exec ../simgs/multiqc.sif multiqc . --verbose --config multiqc_config.yaml 2> multiqc.log

