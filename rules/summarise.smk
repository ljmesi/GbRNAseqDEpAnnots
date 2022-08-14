
#### Join all the data ####
rule join_all_tables:
    input:
        bp = f"{TBLS}/GSEA/genes-GO_BP_Fisher.tsv",
        cc = f"{TBLS}/GSEA/genes-GO_CC_Fisher.tsv",
        mf = f"{TBLS}/GSEA/genes-GO_MF_Fisher.tsv",
        sg_gff = f"{PROC}/annotations/shrinked_padj-filtered_annotated.tsv"
    output:
        report(f"{TBLS}/summarise/combined-genes-GO_Fisher.tsv",
                caption = f"{REP}/summarise/combined-genes-GO_Fisher.rst",
                category = SUM)
    benchmark:
        BMARKS/"summarise"/"join_all_tables.tsv"
    conda:
        f"{ENVS}/summarise.yml"
    log: 
        LOGS/"summarise"/"join_all_table.log"
    script:
        str(SRC/"tbls"/"summarise"/"join_all_tables.R")


rule one_row_per_geneID:
    input:
        f"{TBLS}/summarise/combined-genes-GO_Fisher.tsv"
    output:
        report(f"{TBLS}/summarise/each_gene_in_one_row.tsv",
                caption = f"{REP}/summarise/each_gene_in_one_row.rst",
                category = SUM)
    params:
        str(SRC/"tbls"/"summarise"/"each_gene_in_one_row.py")
    benchmark:
        f"{BMARKS}/summarise/one_row_per_geneID.tsv"
    log:
        f"{LOGS}/summarise/one_row_per_geneID.log"
    conda:
        f"{ENVS}/annotations.yml"
    shell:
        r"""
        sed -i -e 's/GO.ID/GO_ID/g' \
        -e 's/Rank in classic/Rank_in_classic/g' \
        {input}
        python {params} {input} {output} &> {log}
        """

rule interactive_volcano_plot:
    input:
        annotations = f"{TBLS}/summarise/each_gene_in_one_row.tsv",
        de = f"{TBLS}/DE/shrinked_not-filtered.tsv", 
    output:
        fig = f'{FIGS}/summarise/interactive_volcano.html'
    benchmark:
        f"{BMARKS}/summarise/interactive_volcano.tsv"
    log:
        f"{LOGS}/summarise/interactive_volcano.log"
    params:
        str(SRC/"vis"/"summarise"/"interactive_volcano.py")
    conda:
        f"{ENVS}/annotations.yml"
    shell:
        r"""
        python {params} {input.annotations} {input.de} {output} &> {log}
        """
