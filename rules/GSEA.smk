#### Run TopGO and produce output tables and plots ####

rule assign_GO_terms_to_geneIDs:
    input:
        sigGenelist = f"{TBLS}/DE/shrunken_padj-filtered.tsv",
        geneGOmappings = f"{PROC}/annotations/geneIDs_GOs.tsv"
    output:
        f"{TBLS}/GSEA/GO-genes_list_{{GO_category}}_Fisher.tsv"
    params:
        nodeSize = 5
    wildcard_constraints:
        GO_category = "(BP)|(MF)|(CC)"
    conda:
        f"{ENVS}/GSEA.yml"
    benchmark:
        BMARKS/"GSEA"/"assign_GO_terms_to_geneIDs_{{GO_category}}.tsv"
    log: 
        LOGS/"GSEA"/"assign_GO_terms_to_geneIDs_{{GO_category}}.log"
    script:
        str(SRC/"data"/"GSEA"/"assign_GO_terms_to_geneIDs.R")


rule expand_genes:
    input:
        f"{TBLS}/GSEA/GO-genes_list_{{GO_category}}_Fisher.tsv"
    output:
        f"{TBLS}/GSEA/genes-GO_{{GO_category}}_Fisher.tsv"
    wildcard_constraints:
        GO_category = "(BP)|(MF)|(CC)"
    conda:
        f"{ENVS}/GSEA.yml"
    benchmark:
        BMARKS/"GSEA"/"expand_genes_{{GO_category}}.tsv"
    params:
        str(SRC/"data"/"GSEA"/"expandGenes.awk")
    log: 
        LOGS/"GSEA"/"genes_GO_{{GO_category}}.log"
    shell:    
        r"""
        awk \
        -F"\t" \
        -v GOcat={wildcards.GO_category} \
        -f {params} \
        {input} > {output} 2> {log}
        """

rule plot_GO_terms:
    input:
        sigGenelist = f"{TBLS}/DE/shrunken_padj-filtered.tsv",
        geneGOmappings = f"{PROC}/annotations/geneIDs_GOs.tsv"
    output:
        GO_graph = report(f"{FIGS}/GSEA/{{GO_category}}_Fisher_{{algo}}.svg",
                        caption = f"{REP}/GSEA/plot_GO_terms.rst",
                        category = GSEA)
    params:
        nodeSize = 5,
        firstSigNodes = 5
    wildcard_constraints:
        GO_category = "(BP)|(MF)|(CC)"
    conda:
        f"{ENVS}/GSEA.yml"
    benchmark:
        BMARKS/"GSEA"/"{{GO_category}}_plot_GO_terms_{{algo}}.tsv"
    log: 
        LOGS/"GSEA"/"{{GO_category}}_Fisher_{{algo}}.log"
    script:
        str(SRC/"vis"/"GSEA"/"plot_GO_terms.R")
