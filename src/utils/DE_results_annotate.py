#!/usr/bin/env python3
"""
Read DE results tsv file and add annotations from annotations sqlite database
"""

import click
import pandas as pd
from pathlib import Path
import gffutils

def read_gff_db(db_file):
    """Read annotation sqlite database and return FeatureDB instance"""
    click.echo(f"Reading {db_file} database...")
    return gffutils.FeatureDB(db_file)

def match_mRNAs_and_CDSs(mRNA_list, CDS_list):
    """Join together into a tuple all CDS objects in list and their parent mRNA"""
    matched = []
    for mRNA in mRNA_list:
        # Start with an empty list of CDSs for each mRNA
        all_CDSs_of_mRNA = []
        for cds in CDS_list:
            if(mRNA.id == cds.attributes.get("Parent")[0]):
                all_CDSs_of_mRNA.append(cds)
        matched.append((mRNA, all_CDSs_of_mRNA))
    return matched

def collect_all_mRNAs_and_CDSs(db, gene):
    """
    Gather all mRNA and CDS features of a gene into lists and return 
    them as a tuple
    """
    mRNAs = []
    CDSs = []
    for feature in db.children(gene):
        if(feature.featuretype == "mRNA"):
            mRNAs.append(feature)
            continue
        if(feature.featuretype == "CDS"):
            CDSs.append(feature)
            continue
    return (mRNAs, CDSs)

def remove_dublicates(list_of_products):
    """Remove dublicate gene product names from a list"""
    return list(set(list_of_products))

def get_products_from_cds_list(cds_list):
    """Get a list of unique gene product strings from a list of CDS objects"""
    all_products = []
    found_products = False
    for cds in cds_list:
        # If there is a product it is returned inside of a list
        product = cds.attributes.get("Product")
        if product:
            all_products.append(product[0])
            found_products = True
    if found_products:
        return remove_dublicates(all_products)
    else:
        return list()
        
def get_IDs_from_cds_list(cds_list):
    """Get a list of CDS IDs from a list of CDS objects"""
    return [cds.id for cds in cds_list]

def create_gene_annotations(db)
    """Create gene annotions dictionary with all essential info in memory"""
    click.echo("Creation of annotions dictionary started...")
    gene_annotations = {}
    for gene in db.features_of_type('gene'):
        # Gather the gene's all mRNAs and CDSs into own lists
        mRNAs, CDSs = collect_all_mRNAs_and_CDSs(db, gene)
        # Iterate over CDSs and their parent mRNA and store all relevant info into
        # gene_annotation dictionary
        for tx_cds in match_mRNAs_and_CDSs(mRNAs, CDSs):        
            mRNA = tx_cds[0]
            all_CDSs_of_mRNA = tx_cds[1]
            cds_products = get_products_from_cds_list(all_CDSs_of_mRNA)
            cds_IDs = get_IDs_from_cds_list(all_CDSs_of_mRNA)
            # Populate gene_annotations dictionary
            if gene.id in gene_annotations:
                gene_annotations[gene.id].append({
                    "Gene_name" : gene.attributes.get("Name"),
                    "Transcript_ID" : mRNA.id,
                    "Ontology_term" : mRNA.attributes.get("Ontology_term"),
                    "Dbxref" : mRNA.attributes.get("Dbxref"),
                    "CDS_ID" : cds_IDs,
                    "CDS_Product" : cds_products
                    })
            else:
                gene_annotations[gene.id] = [{
                    "Gene_name" : gene.attributes.get("Name"),
                    "Transcript_ID" : mRNA.id,
                    "Ontology_term" : mRNA.attributes.get("Ontology_term"),
                    "Dbxref" : mRNA.attributes.get("Dbxref"),
                    "CDS_ID" : cds_IDs,
                    "CDS_Product" : cds_products
                    }]
    click.echo("Gene annotions dictionary finished")
    return gene_annotations

@click.command()
@click.argument('db-file', type=click.Path(exists=True))
@click.argument('unannotated-DE-results', type=click.Path(exists=True))
@click.argument('annotated-DE-results', type=click.Path())
def main(db_file, unannotated_DE_results, annotated_DE_results):
    """Run the main CLI application"""
    db = read_gff_db(db_file)
    click.echo(f"{db_file} database read")
    gene_annotations = create_gene_annotations(db)

if __name__ == '__main__':
    main()