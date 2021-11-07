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
        expand(f"{FIGS}/DE/normalised.{{eda_type}}.{{transformation_type}}.svg",
               eda_type = ["correlation","pca"],
               transformation_type = ["rlog","vs"]),
        # Tables containing DE results
        expand(f"{TBLS}/DE/{{shrink_status}}_{{subsetting}}.tsv",
                shrink_status = SHRINKAGE,
                subsetting = ["not-filtered","padj-filtered"]),

        # DESeq Model QC plots
        expand(f"{FIGS}/DE/{{qc_type}}.svg",
                qc_type = ["log_mean_log_variance","deseq_dispersions"]),
        # Outliers recognised by DESeq2 gather in one table
        f"{TBLS}/DE/outliers.tsv",
        # Visualising the DE results
        expand(f"{FIGS}/DE/MAplot_{{shrink_status}}.svg",
                shrink_status = SHRINKAGE),
        f"{FIGS}/DE/volcanoplot.svg",
        f"{FIGS}/DE/DE_heatmap.svg",
        f"{PROC}/DE/genelist_padj-filtered.tsv"
        

        

                
               
#### Analysis modules ####

include: "rules/DE.smk"
