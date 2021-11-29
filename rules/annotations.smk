
#### Preprocess annotations useful in the next steps ####

rule create_annotations_db:
    input:
        f"{REF}/genome.genes.gff3"
    output:
        protected(f"{REF}/genome.genes.sqlite")
    params:
        script = str(SRC/"tbls"/"annotations"/"create_annotations_db.py"),
        logging = config["software"]["utils"]["csvkit_log"]
    conda:
        f"{ENVS}/annotations.yml"
    benchmark:
        BMARKS/"annotations"/"create_annotations_db.tsv"
    log:
        LOGS/"annotations"/"create_annotations_db.log"
    shell:
        r"""
        python {params.script} {input} {output} &> {log}
        """

rule create_gene_universe:
    input:
        f"{REF}/genome.genes.sqlite"
    output:
        f"{PROC}/annotations/geneIDs_GOs.tsv"
    params:
        str(SRC/"tbls"/"annotations"/"create_gene_universe.py")
    conda:
        f"{ENVS}/annotations.yml"
    benchmark:
        BMARKS/"annotations"/"create_gene_universe.tsv"
    log:
        LOGS/"annotations"/"create_gene_universe.log"
    shell:
        r"""
        python {params} {input} {output} &> {log}
        """

#### Add annotations from gff file to DE genes ####

rule add_gff_annotations:
    input:
        sqlite = f"{REF}/genome.genes.sqlite",
        de_data = f"{TBLS}/DE/shrinked_padj-filtered.tsv"
    output:
        f"{PROC}/annotations/gff_DE_geneIDs.tsv"
    params:
        str(SRC/"tbls"/"annotations"/"add_gff_annotations.py")
    conda:
        f"{ENVS}/annotations.yml"
    benchmark:
        BMARKS/"annotations"/"add_gff_annotations.tsv"
    log:
        LOGS/"annotations"/"add_gff_annotations.log"
    shell:
        r"""
        python {params} {input.sqlite} {input.de_data} {output} &> {log}
        """