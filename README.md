# *Gerris buenoi* RNAseq Differential Expression analysis + Annotations
- [*Gerris buenoi* RNAseq Differential Expression analysis + Annotations](#gerris-buenoi-rnaseq-differential-expression-analysis--annotations)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Gene Set Enrichment Analysis](#gene-set-enrichment-analysis)
  - [Summarisation](#summarisation)
  - [Makefile](#makefile)
  - [Raw, Intermediary and Reference data](#raw-intermediary-and-reference-data)
  - [Results](#results)

The script used for indexing *Gerris buenoi* reference genome using [STAR v. 2.7.9a](https://github.com/alexdobin/STAR) can be found in [`exp/submit_star_alignments/STARindex.sh`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/exp/submit_star_alignments/STARindex.sh). 
The scripts for performing mapping of reads onto the reference genome can be found in: 
- [`exp/submit_star_alignments/STARmap.sh`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/exp/submit_star_alignments/STARmap.sh)
- [`exp/submit_star_alignments/submit_star_alignmets.sh`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/exp/submit_star_alignments/submit_star_alignmets.sh)
The latter script uses [GNU parallel](https://www.gnu.org/software/parallel/) [v. 20180822](https://doi.org/10.5281/zenodo.1146014) which takes as input a simple list of sample names found in [`data/ref/sample_names.txt`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/data/ref/sample_names.txt).

## Differential Expression Analysis

The script for running differential expression (DE) analysis with [DESeq2 v. 1.3.4](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) can be found in [`src/tbls/DE/execute_DE.R`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/src/tbls/DE/execute_DE.R) and a [snakemake v. 7.12.1](https://snakemake.readthedocs.io/en/stable/) script for running the R script can be found in [`rules/DE.smk`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/rules/DE.smk). 

## Gene Set Enrichment Analysis

The script for producing gene universe from genome annotation used in gene set enrichment analysis (GSEA) performed with [topGO R package](https://bioconductor.org/packages/release/bioc/html/topGO.html) v. 2.46.0 can be found in [`src/tbls/annotations/create_gene_universe.py`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/src/tbls/annotations/create_gene_universe.py) and was run with [`rules/annotations.smk`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/rules/annotations.smk). GSEA analysis was executed with [`src/data/GSEA/assign_GO_terms_to_geneIDs.R`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/src/data/GSEA/assign_GO_terms_to_geneIDs.R) and run with [`rules/GSEA.smk`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/rules/GSEA.smk).

## Summarisation

Finally, significant (<img src="https://render.githubusercontent.com/render/math?math=\text{adjusted p-value} \leq 0.01">) results from DE analysis, GSEA and genome annotations were combined into one table and interactive [plotly (v.5.9.0)](https://plotly.com/python/) volcano plot. These were executed with [`src/tbls/summarise/each_gene_in_one_row.py`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/src/tbls/summarise/each_gene_in_one_row.py) and [`src/vis/summarise/interactive_volcano.py`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/src/vis/summarise/interactive_volcano.py) respectively.

## Makefile

The root of the repository contains a [`Makefile`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/Makefile) which provides shortcut commands for the most common tasks that can be executed in the project. The make rules are set to be run inside a conda environment defined by [`envs/snakemake.yml`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/envs/snakemake.yml).

## Raw, Intermediary and Reference data

The [`data`](https://github.com/ljmesi/GbRNAseqDEpAnnots/tree/main/data) directory contains intermediary, raw and reference files used or produced by the snakemake pipeline. The genome annotation in [`data/ref`](https://github.com/ljmesi/GbRNAseqDEpAnnots/tree/main/data/ref) was gitignored and thus is not included in the repository. A table with read counts mapping to genes used as input for differential expression analysis can be found in [`data/raw/all_ReadsPerGene.out.tsv`](https://github.com/ljmesi/GbRNAseqDEpAnnots/blob/main/data/raw/all_ReadsPerGene.out.tsv).

## Results

All the results for the analyses can be found in [`results`](https://github.com/ljmesi/GbRNAseqDEpAnnots/tree/main/results) directory. 
