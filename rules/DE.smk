
#### Load samples and sample QC ####
rule read_normalise_counts:
    input:
        raw_counts_table = f"{RAW}/STAR_counts_all.csv",
        experiment_metadata = config["DE"]["Experimental_design"]
    output:
        DESeq_obj = f"{PROC}/DE/DESeq_raw.RDS",
        raw_counts = f"{PROC}/DE/raw_counts.tsv",
        normalised_counts_DESeq = f"{PROC}/DE/DESeq_normalised.RDS",
        normalised_counts = f"{PROC}/DE/normalised_counts.tsv"
    conda:
        f"{ENVS}DE.yml"
    benchmark:
        f"{BMARKS}/DE/read_normalise_counts.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        log = f"{LOGS}/DE/read_normalise_counts.log"
    script:
        f"{SRC}/data/DE/read_normalise_counts.R"

