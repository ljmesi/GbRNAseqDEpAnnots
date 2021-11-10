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
               transformation_type = ["rlog","vs"])
               
#### Analysis modules ####

include: "rules/DE.smk"
