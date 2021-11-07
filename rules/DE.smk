
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
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"read_normalise_counts.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"read_normalise_counts.log"
    script:
        str(SRC/"data"/"DE"/"read_normalise_counts.R")


# RNA-seq counts sample-wise QC and EDA
rule plot_correlation_heatmaps_from_transf_DESeq:
    input:
        normalised_counts_DESeq = f"{PROC}/DE/DESeq_normalised.RDS",
        metadata = config["DE"]["Experimental_design"]
    output:
        rlog_transformed_corr_heatmap = report(f"{FIGS}/DE/normalised.correlation.rlog.svg",
                                               caption = f"{REP}/DE/rlog_correlation_heatmap.rst",
                                               category = QC_EDA),
        vs_transformed_corr_heatmap = report(f"{FIGS}/DE/normalised.correlation.vs.svg",
                                             caption = f"{REP}/DE/vs_correlation_heatmap.rst",
                                             category = QC_EDA)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_correlation_heatmaps_from_transf_DESeq.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"transformed.correlation.heatmap.log"
    script:
        str(SRC/"vis"/"DE"/"plot_correlation_heatmap_transform_DESeq.R")

rule plot_pca_from_transf_DESeq:
    input:
        normalised_counts_DESeq = f"{PROC}/DE/DESeq_normalised.RDS",
        metadata = config["DE"]["Experimental_design"]
    output:
        rlog_transformed_pca = report(f"{FIGS}/DE/normalised.pca.rlog.svg",
                                               caption = f"{REP}/DE/rlog_qc_pca.rst",
                                               category = QC_EDA),
        vs_transformed_pca = report(f"{FIGS}/DE/normalised.pca.vs.svg",
                                             caption = f"{REP}/DE/vs_qc_pca.rst",
                                             category = QC_EDA)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_pca_from_transf_DESeq.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"qc.pca.log"
    script:
        str(SRC/"vis"/"DE"/"plot_pca_transform_DESeq.R")


#### DE analysis ####

rule execute_DE:
    input:
        DESeq_raw = f"{PROC}/DE/DESeq_raw.RDS"
    output:
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS",
        DESeq_results_non_shrinked = f"{PROC}/DE/not-shrinked.RDS",
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
        tbl_non_shrinked = f"{TBLS}/DE/not-shrinked_not-filtered.tsv",
        tbl_non_shrinked_padj_filtered = f"{TBLS}/DE/not-shrinked_padj-filtered.tsv",
        tbl_shrinked = report(f"{TBLS}/DE/shrinked_not-filtered.tsv", 
                                caption = f"{REP}/DE/shrinked.rst",
                                category = DE),
        tbl_shrinked_padj_filtered = report(f"{TBLS}/DE/shrinked_padj-filtered.tsv", 
                                    caption = f"{REP}/DE/shrinked_sig.rst",
                                    category = DE),
    params:
        p_adj_limit = padj_limit,
        contrast = config["DE"]["contrasts"],
        pAdjustMethod = config["DE"]["pAdjustMethod"]
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"execute_DE.tsv"
    threads:
        config["software"]["threads"]["DEseq2"]
    log:
        LOGS/"DE"/"execute_DE.log"
    script:
        str(SRC/"tbls"/"DE"/"execute_DE.R")

