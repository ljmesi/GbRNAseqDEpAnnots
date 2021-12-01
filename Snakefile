from pathlib import Path
from snakemake.utils import min_version, validate
import pandas as pd
##### set minimum snakemake version #####
min_version("6.10.0")

configfile: "conf/config.yml"
#### Load report description #### 
report: "report/workflow.rst"

include: "rules/common.smk"

rule all:
    input:
        #### Differential expression with DESeq2 ####
        # Sample QC and EDA
        f"{PROC}/DE/normalised_counts.tsv",
        expand(f"{FIGS}/DE/{{eda_type}}.{{transformation_type}}.svg",
                eda_type = ["gene_clustering","pca"],
                transformation_type = ["rlog","vs"]),
        # Tables containing DE results
        expand(f"{TBLS}/DE/{{shrink_status}}_{{subsetting}}.tsv",
                shrink_status = SHRINKAGE,
                subsetting = ["not-filtered","padj-filtered"]),
        expand(f"{FIGS}/DE/{{transformation}}_{{plot}}.svg",
                transformation = ["rlog","vs"],
                plot = ["count_matrix","sample_to_sample_distances"]),

        # DESeq Model QC plots
        f"{FIGS}/DE/log10_cooks_distances.svg",
        f"{FIGS}/DE/deseq_dispersions.svg",
        # Outliers recognised by DESeq2 gather in one table
        f"{TBLS}/DE/outliers.tsv",
        # Visualising the DE results
        expand(f"{FIGS}/DE/MAplot_{{shrink_status}}.svg",
                shrink_status = SHRINKAGE),
        f"{FIGS}/DE/volcanoplot.svg",
        f"{FIGS}/DE/DE_heatmap.svg",

        #### gff-annotations ####
        f"{REF}/genome.genes.sqlite",
        f"{PROC}/annotations/geneIDs_GOs.tsv",
        f"{PROC}/annotations/gff_DE_geneIDs.tsv",

        #### GSEA of DE results with topGO ####
        expand(f"{TBLS}/GSEA/GO-genes_list_{{GO_category}}_Fisher.tsv",
               GO_category = config["GSEA"]["GOcategories"]),
        expand(f"{TBLS}/GSEA/genes-GO_{{GO_category}}_Fisher.tsv",
               GO_category = config["GSEA"]["GOcategories"]),
        expand(f"{FIGS}/GSEA/{{GO_category}}_Fisher_{{algo}}.svg",
               GO_category = config["GSEA"]["GOcategories"],
               algo = config["GSEA"]["GOhierarchyTraversingAlgorithms"]),
        # f"{TBLS}/summarise/combined-genes-GO_Fisher.tsv",

        #### Final parts of pipeline (join all data into one large table) ####
        f"{PROC}/summarise/combined-genes-GO_Fisher.tsv",
        f"{TBLS}/summarise/all_compressed.tsv"

#### Analysis modules ####

include: "rules/DE.smk"
include: "rules/annotations.smk"
include: "rules/GSEA.smk"
include: "rules/summarise.smk"
