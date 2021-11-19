
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


#RNA-seq counts sample-wise QC and EDA

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
                                            category = QC_EDA),
        rlog_transformed = f"{PROC}/DE/rlog_transformed.RDS",
        vs_transformed = f"{PROC}/DE/vs_transformed.RDS"
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
        rlog_transformed = f"{PROC}/DE/rlog_transformed.RDS",
        vs_transformed = f"{PROC}/DE/vs_transformed.RDS",
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


#### DESeq Model QC ####

rule plot_DESeq_log_mean_log_variance:
    input:
        raw_counts = f"{PROC}/DE/raw_counts.tsv"
    output:
        pl = report(f"{FIGS}/DE/log_mean_log_variance.svg",
                    caption = f"{REP}/DE/log_mean_log_variance.rst",
                    category = QC_DESEQ)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_DESeq_log_mean_log_variance.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"plot_log_mean_log_variance.log"
    script:
        str(SRC/"vis"/"DE"/"plot_log_mean_log_variance.R")

rule plot_dispersions:
    input:
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS"
    output:
        dispersion_plot = report(f"{FIGS}/DE/deseq_dispersions.svg",
                                    caption = f"{REP}/DE/dispersions.rst",
                                    category = QC_DESEQ)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_dispersions.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"plot_dispersions.log"
    script:
        str(SRC/"vis"/"DE"/"plot_dispersions.R")

rule find_outliers:
    input:
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS",
        raw_counts = f"{PROC}/DE/raw_counts.tsv",
        metadata = config["DE"]["Experimental_design"]
    output:
        table_outliers = report(f"{TBLS}/DE/outliers.tsv", 
                                caption = f"{REP}/DE/outliers.rst",
                                category = DE)
    params:
        p_adj_limit = padj_limit,
        contrast = config["DE"]["contrasts"],
        pAdjustMethod = config["DE"]["pAdjustMethod"]
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"find_outliers.tsv"
    threads: 
        config["software"]["threads"]["DEseq2"]
    log: 
        LOGS/"DE"/"outliers.log"
    script:
        str(SRC/"tbls"/"DE"/"find_outliers.R")

#### Visualisation of differential expression ####

rule plot_MAplot:
    input:
        DESeq_results = f"{PROC}/DE/{{shrink_status}}.RDS"
    output:
        MAplot = report(f"{FIGS}/DE/MAplot_{{shrink_status}}.svg",
                        caption = f"{REP}/DE/MAplots.rst",
                        category = DE)
    wildcard_constraints:
        shrink_status="(not-shrinked)|(shrinked)"
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_MAplot_{{shrink_status}}.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"MAplot_{{shrink_status}}.log"
    script:
        str(SRC/"vis"/"DE"/"plot_MAplot.R")


rule plot_volcano:
    input:
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS"
    output:
        volcano_plot = report(f"{FIGS}/DE/volcanoplot.svg",
                                caption = f"{REP}/DE/volcanoplot.rst",
                                category = DE)
    params:
        p_adj_limit = padj_limit,
        volcano_top_genes = config["DE"]["volcano_top_genes"]
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_volcano.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"plot_volcano.log"
    script:
        str(SRC/"vis"/"DE"/"plot_volcano.R")


rule plot_DE_heatmap:
    input:
        normalised_counts = f"{PROC}/DE/normalised_counts.tsv",
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
        metadata = config["DE"]["Experimental_design"]
    output:
        DE_heatmap = report(f"{FIGS}/DE/DE_heatmap.svg",
                            caption = f"{REP}/DE/DE_heatmap.rst",
                            category = DE)
    params:
        p_adj_limit = padj_limit
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_DE_heatmap.tsv"
    threads: 
        config["software"]["threads"]["default"]
    log: 
        LOGS/"DE"/"DE_heatmap.log"
    script:
        str(SRC/"vis"/"DE"/"plot_DE_heatmap.R")

#### Prepare a list of genes of interest ####

rule select_geneIDs_and_padj:
    input:
        f"{TBLS}/DE/shrinked_padj-filtered.tsv"
    output:
        f"{PROC}/DE/genelist_padj-filtered.tsv"
    params:
        logging = config["software"]["utils"]["csvkit_log"]
    threads: 
        config["software"]["threads"]["default"]
    conda:
        f"{ENVS}/GSEA.yml"
    benchmark:
        BMARKS/"DE"/"select_geneIDs_and_padj.tsv"
    log: 
        LOGS/"DE"/"select_geneIDs_and_padj.log"
    shell:
        r"""
        cut -f1,7 {input} > {output} 2> {log} &&
        bash {params.logging} "-H" {output} &>> {log}
        """

