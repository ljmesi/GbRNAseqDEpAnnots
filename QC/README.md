# Directory contents

This directory contains quality control related files.

## FastQC

FastQC was run with commands:

```shell
cd /proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/QC
mkdir fastqc
module load bioinfo-tools gnuparallel/20180822
cat ../data/ref/sample_names.txt | \
parallel --eta \
--progress \
--bar \
--verbose \
--jobs 8 \
"singularity exec ../simgs/fastqc.sif fastqc ../data/raw/{}_R1_001.fastq.gz ../data/raw/{}_R2_001.fastq.gz --format fastq --outdir fastqc"
```

## Seqkit stats

Seqkit stats was run with commands:

```shell
module load bioinfo-tools gnuparallel/20180822
cat ../data/ref/sample_names.txt | \
parallel --eta \
--keep-order \
--interactive \
--progress \
--bar \
--verbose \
--jobs 8 \
-N1 \
"singularity exec ../simgs/seqkit.sif seqkit stats ../data/raw/{}_R1_001.fastq.gz ../data/raw/{}_R2_001.fastq.gz >> seqkit_stats.tsv 2>&1"
```

## MultiQC

In order to summarise the initial FastQC results MultiQC was run:

```shell
singularity exec ../simgs/multiqc.sif \
multiqc . \
--verbose 2> multiqc.log
```
