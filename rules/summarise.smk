
#### Join all the data ####
rule join_all_tables:
    input:
        BP = f"{TBLS}/GSEA/genes-GO_BP_Fisher.tsv",
        CC = f"{TBLS}/GSEA/genes-GO_CC_Fisher.tsv",
        MF = f"{TBLS}/GSEA/genes-GO_MF_Fisher.tsv",
        SG = f"{TBLS}/DE/shrinked_padj-filtered.tsv",
        GFF = f"{PROC}/annotations/gff_DE_geneIDs.tsv"
    output:
        f"{PROC}/summarise/combined-genes-GO_Fisher.tsv"
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
        f"{PROC}/summarise/combined-genes-GO_Fisher.tsv"
    output:
        report(f"{TBLS}/summarise/all_compressed.tsv",
               caption = f"{REP}/summarise/tsv_compressed.rst",
               category = SUM)
    params:
        str(SRC/"tbls"/"summarise"/"compress_to_one_row.py")    
    benchmark:
        f"{BMARKS}/summarise/one_row_per_geneID.tsv"
    log:
        f"{LOGS}/summarise/one_row_per_geneID.log"
    conda:
        f"{ENVS}/annotations.yml"
    shell:
        r"""
        sed -i -e 's/GO.ID/GO_ID/g' \
        -e 's/Rank in classicFisher/Rank_in_classicFisher/g' \
        {input}
        python {params} {input} {output} &> {log}
        """
