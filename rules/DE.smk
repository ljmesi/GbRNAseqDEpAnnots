
#### Load samples and sample QC ####
rule read_normalise_counts:
    input:
        raw_counts_table = f"{RAW}/all_ReadsPerGene.out.tsv",
        experiment_metadata = config["DE"]["Experimental_design"]
    output:
        DESeq_obj = f"{PROC}/DE/DESeq_raw.RDS",
        normalised_counts = f"{PROC}/DE/normalised_counts.tsv"
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"read_normalise_counts.tsv"
    log: 
        LOGS/"DE"/"read_normalise_counts.log"
    script:
        str(SRC/"data"/"DE"/"read_normalise_counts.R")


#### QC and EDA ####

rule plot_gene_clustering:
    input:
        DESeq_counts = f"{PROC}/DE/DESeq_raw.RDS"
    output:
        rlog_gene_clustering_heatmap = report(f"{FIGS}/DE/gene_clustering.rlog.svg",
                                                caption = f"{REP}/DE/rlog_gene_clustering_heatmap.rst",
                                                category = QC),
        vs_gene_clustering_heatmap = report(f"{FIGS}/DE/gene_clustering.vs.svg",
                                            caption = f"{REP}/DE/vs_gene_clustering_heatmap.rst",
                                            category = QC),
        # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#effects-of-transformations-on-the-variance
        rlog_meanSdPlot = f"{FIGS}/DE/rlog_meanSdPlot.svg",
        vs_meanSdPlot = f"{FIGS}/DE/vs_meanSdPlot.svg",
        # Transformed data for further use
        rlog_transformed = f"{PROC}/DE/rlog_transformed.RDS",
        vs_transformed = f"{PROC}/DE/vs_transformed.RDS"
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_gene_clustering.tsv"
    log: 
        LOGS/"DE"/"plot_gene_clustering.log"
    script:
        str(SRC/"vis"/"DE"/"plot_gene_clustering_transform_DESeq.R")


rule plot_vs_QC_heatmap:
    input:
        f"{PROC}/DE/vs_transformed.RDS"
    output:
        f"{FIGS}/DE/vs_qc_sample_to_sample_distances.svg"
    params:
        num_highest_means = 50
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_vs_QC_heatmap.tsv"
    log: 
        LOGS/"DE"/"plot_vs_QC_heatmap.log"
    script:
        str(SRC/"vis"/"DE"/"plot_vs_QC_heatmap.R")


rule plot_QC_heatmaps:
    input:
        f"{PROC}/DE/{{transformation}}_transformed.RDS"
    output:
        # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
        count_matrix = report(f"{FIGS}/DE/{{transformation}}_count_matrix.svg", 
                                caption = f"{REP}/DE/count_matrix.rst", 
                                category = QC),
        # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances
        sample_to_sample_distances = report(f"{FIGS}/DE/{{transformation}}_sample_to_sample_distances.svg", 
                                            caption = f"{REP}/DE/sample_to_sample_distances.rst", 
                                            category = QC)
    params:
        num_highest_means = 50
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"{transformation}_plot_QC_heatmaps.tsv"
    log: 
        LOGS/"DE"/"{transformation}_plot_QC_heatmaps.log"
    script:
        str(SRC/"vis"/"DE"/"plot_QC_heatmaps.R")


rule plot_pca_from_transf_DESeq:
    input:
        rlog_transformed = f"{PROC}/DE/rlog_transformed.RDS",
        vs_transformed = f"{PROC}/DE/vs_transformed.RDS",
        metadata = config["DE"]["Experimental_design"]
    output:
        rlog_transformed_pca = report(f"{FIGS}/DE/pca.rlog.svg",
                                        caption = f"{REP}/DE/rlog_qc_pca.rst",
                                        category = QC),
        vs_transformed_pca = report(f"{FIGS}/DE/pca.vs.svg",
                                    caption = f"{REP}/DE/vs_qc_pca.rst",
                                    category = QC)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_pca_from_transf_DESeq.tsv"
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
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
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
        config["threads"]["DEseq2"]
    log:
        LOGS/"DE"/"execute_DE.log"
    script:
        str(SRC/"tbls"/"DE"/"execute_DE.R")


#### DESeq Model QC ####

rule plot_log10_cooks_distances:
    input:
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS"
    output:
        report(f"{FIGS}/DE/log10_cooks_distances.svg",
                caption = f"{REP}/DE/log10_cooks_distances.rst",
                category = QC)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_log10_cooks_distances.tsv"
    log: 
        LOGS/"DE"/"plot_log10_cooks_distances.log"
    script:
        str(SRC/"vis"/"DE"/"plot_log10_cooks_distances.R")

rule plot_dispersions:
    input:
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS"
    output:
        report(f"{FIGS}/DE/deseq_dispersions.svg",
                caption = f"{REP}/DE/dispersions.rst",
                category = QC)
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_dispersions.tsv"
    log: 
        LOGS/"DE"/"plot_dispersions.log"
    script:
        str(SRC/"vis"/"DE"/"plot_dispersions.R")

rule find_outliers:
    input:
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
        DESeq_analysis_obj = f"{PROC}/DE/DESeq-obj.RDS",
        raw_counts = f"{RAW}/all_ReadsPerGene.out.tsv",
        metadata = config["DE"]["Experimental_design"]
    output:
        table_outliers = report(f"{TBLS}/DE/outliers.tsv", 
                                caption = f"{REP}/DE/outliers.rst",
                                category = QC)
    params:
        p_adj_limit = padj_limit,
        contrast = config["DE"]["contrasts"],
        pAdjustMethod = config["DE"]["pAdjustMethod"]
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"find_outliers.tsv"
    threads: 
        config["threads"]["DEseq2"]
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
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_MAplot_{shrink_status}.tsv"
    log: 
        LOGS/"DE"/"MAplot_{shrink_status}.log"
    script:
        str(SRC/"vis"/"DE"/"plot_MAplot.R")


rule plot_volcano:
    input:
        DESeq_results_shrinked = f"{PROC}/DE/shrinked.RDS",
        gene_annotations = f"{TBLS}/summarise/each_gene_in_one_row.tsv"
    output:
        volcano_plot = f"{FIGS}/DE/volcanoplot_sketch.svg",
        volcano_plot_pub = report(f"{FIGS}/DE/volcanoplot.svg",
                                caption = f"{REP}/DE/volcanoplot.rst",
                                category = DE)
    params:
        p_adj_limit = padj_limit
    conda:
        f"{ENVS}/DE.yml"
    benchmark:
        BMARKS/"DE"/"plot_volcano.tsv"
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
    log: 
        LOGS/"DE"/"DE_heatmap.log"
    script:
        str(SRC/"vis"/"DE"/"plot_DE_heatmap.R")

