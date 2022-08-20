# *Gerris buenoi* RNAseq Differential Expression analysis + Annotations
- [*Gerris buenoi* RNAseq Differential Expression analysis + Annotations](#gerris-buenoi-rnaseq-differential-expression-analysis--annotations)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Gene Set Enrichment Analysis](#gene-set-enrichment-analysis)
  - [Summarisation](#summarisation)
  - [Raw, Intermediary and Reference data](#raw-intermediary-and-reference-data)
  - [Results](#results)

The script used for indexing *Gerris buenoi* reference genome using [STAR v. 2.7.9a](https://github.com/alexdobin/STAR) can be found in `exp/submit_star_alignments/STARindex.sh`. 
The scripts for performing mapping of reads onto the reference genome can be found in: 
- `exp/submit_star_alignments/STARmap.sh`
- `exp/submit_star_alignments/submit_star_alignmets.sh`
The latter script uses [GNU parallel](https://www.gnu.org/software/parallel/) [v. 20180822](https://doi.org/10.5281/zenodo.1146014) which takes as input a simple list of sample names found in `exp/submit_star_alignments/sample_names.txt`.

## Differential Expression Analysis

The script for running differential expression (DE) analysis with [DESeq2 v. 1.3.4](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) can be found in `src/tbls/DE/execute_DE.R` and a [snakemake v. 7.12.1](https://snakemake.readthedocs.io/en/stable/) script for running the R script can be found in `rules/DE.smk`. 

## Gene Set Enrichment Analysis

The script for producing gene universe from genome annotation used in gene set enrichment analysis (GSEA) performed with [topGO R package](https://bioconductor.org/packages/release/bioc/html/topGO.html) v. 2.46.0 can be found in `src/tbls/annotations/create_gene_universe.py` and was run with `rules/annotations.smk`. GSEA analysis was executed with `src/data/GSEA/assign_GO_terms_to_geneIDs.R` and run with `rules/GSEA.smk`.

## Summarisation

Finally, significant ($\text{adjusted p-value} \leq 0.01$) results from DE analysis, GSEA and genome annotations were combined into one table and interactive [plotly (v.5.9.0)](https://plotly.com/python/) volcano plot. These were executed with `src/tbls/summarise/each_gene_in_one_row.py` and `src/vis/summarise/interactive_volcano.py` respectively.

## Raw, Intermediary and Reference data

The `data` directory contains intermediary, raw and reference files used or produced by the snakemake pipeline. The genome annotation in `data/ref` was gitignored and thus is not included in the repository. A table with read counts mapping to genes used as input for differential expression analysis can be found in `data/raw/all_ReadsPerGene.out.tsv`.

## Results

All the results for the analyses can be found in `results` directory. 
