#!/usr/bin/env python3
"""
Read DE results tsv file and add annotations from annotations sqlite database
"""

import click
import pandas as pd
from pathlib import Path
import gffutils
import csv
from collections import namedtuple

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

def create_gene_annotations(db):
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


def parse_refs(ref):
    """Read an unparsed reference element and return db id"""
    return (":").join(ref.split(":")[1:])

def add_ref(ref_dict, reference_db_name, unparsed_ref_element):
    """Create a list of reference values"""
    if reference_db_name in ref_dict:
        ref_dict[reference_db_name].append(parse_refs(unparsed_ref_element))
    else:
        ref_dict[reference_db_name] = [parse_refs(unparsed_ref_element)]
    return ref_dict

def parse_dbxref_list(dbxref_list, all_found):
    """
    Traverse a list of unparsed db references and add them to a dict of lists
    Keep also track of which db references didn't exist in the list
    """
    dbxref = {}
    for unparsed_reference in dbxref_list:
        if (unparsed_reference.startswith("UniProtKB") or unparsed_reference.startswith("Swiss-Prot")):
            dbxref = add_ref(dbxref, "UniProtKB", unparsed_reference)
            all_found['UniProtKB'] = True
        elif (unparsed_reference.startswith("GeneID")):
            dbxref = add_ref(dbxref, "GeneID", unparsed_reference)
            all_found['GeneID'] = True
        elif (unparsed_reference.startswith("KEGG")):
            dbxref = add_ref(dbxref, "KEGG", unparsed_reference)
            all_found['KEGG'] = True
        elif (unparsed_reference.startswith("PFAM")):
            dbxref = add_ref(dbxref, "PFAM", unparsed_reference)
            all_found['PFAM'] = True
        elif (unparsed_reference.startswith("InterPro")):
            dbxref = add_ref(dbxref, "InterPro", unparsed_reference)
            all_found['InterPro'] = True
        elif (unparsed_reference.startswith("EMBL")):
            dbxref = add_ref(dbxref, "EMBL", unparsed_reference)
            all_found['EMBL'] = True
        else:
            # If not matched with any, create an unparsed unknown element
            dbxref["Unknown"] = unparsed_reference
            all_found['Unknown'] = True
            click.echo(f"Unknown reference added: {unparsed_reference}")
    return (dbxref, all_found)

def ntuple_to_list(ntuple):
    """Convert named tuple values to list"""
    return list(ntuple)

def parse_Dbxref(dbxref_list):
    """Parse a list dbxrefs and return a dict of db as key and ids as values"""
    dbxref = {}
    # Handle completely missing Dbxref data
    if dbxref_list is None:
        return {
            "UniProtKB":[""],
            "GeneID":[""],
            "KEGG":[""],
            "PFAM":[""],
            "InterPro":[""],
            "EMBL":[""]
            }
    
    # Keep track of partially missing Dbxref data
    all_found = {
        "UniProtKB" : False,
        "GeneID" : False,
        "KEGG" : False,
        "PFAM" : False,
        "InterPro" : False,
        "EMBL" : False,
        "Unknown" : False
    }
    # Create dicts for db ids and how partially the references were found
    dbxref, all_found = parse_dbxref_list(dbxref_list, all_found)
    
    # Fill in missing data with empty strings in a list
    for db, found in all_found.items():
        if not found:
            dbxref[db] = [""]
    return dbxref

def remove_duplicate_values(dbxref_dict):
    """Go through all db references and remove duplicate db IDs"""
    for key, value in dbxref_dict.items():
        dbxref_dict[key] = remove_dublicates(value)
    return dbxref_dict

def convert_nones_to_list(element):
    """Convert NoneTypes are converted to empty lists"""
    if element is None:
        return []
    else:
        return element


def concat_list(in_list,concat_char = ","):
    """Concatenate all list elements and return a single element list of the results"""
    if isinstance(in_list, list):
        return [concat_char.join(in_list)]
    else:
        return [""]
    

@click.command()
@click.argument('db-file', type=click.Path(exists=True))
@click.argument('unannotated-de', type=click.Path(exists=True))
@click.argument('annotated-de', type=click.Path())
def main(db_file, unannotated_de, annotated_de):
    """Run the main CLI application"""
    db = read_gff_db(db_file)
    click.echo(f"{db_file} database read")
    gene_annotations = create_gene_annotations(db)
    
    # Read and write gff annotations
    with open(unannotated_de, newline="") as infile, open(annotated_de, 'wt') as out_file:
        # Change path objects to strings
        unannotated_fn = str(unannotated_de)
        annotated_fn = str(annotated_de)

        # Initialise reading and writing tsvs
        reader = csv.reader(infile, delimiter="\t")
        click.echo(f"Opened {unannotated_fn} file for reading")
        tsv_writer = csv.writer(out_file, delimiter='\t')
        click.echo(f"Opened {annotated_fn} file for writing")

        # Write tsv header
        gff_columns = ["Gene_name", "GeneID", "CDS_Products", "Gff_ontology_terms", "Transcript_IDs", "UniProtKB_Swiss_Prot", "KEGG", "PFAM", "InterPro", "EMBL"]
        tsv_writer.writerow(["Genes"]+gff_columns)
        click.echo(f"Header written to {annotated_fn}")

        # Read tsv file row-wise
        Data = namedtuple("Data", next(reader))  # get names from column headers
        click.echo(f"Starting to write the parsed gff-data to {annotated_fn} file")
        for data in map(Data._make, reader):
            # Capture DE data into a list
            data_row = ntuple_to_list(data)
            # Reset all lists/dict for each gene
            ontology_terms = []
            transcript_IDs = []
            CDS_Products = []
            dbxref_dict = {}
            
            # Gather all annotation data 
            # Gene annotations for each gene ID consist of a list of one or more 
            # transcripts with all available annotations
            for transcript in gene_annotations[data.Genes]:
                # gene_name is the same for all the transcripts
                gene_name = transcript["Gene_name"]
                transcript_IDs.append(transcript["Transcript_ID"])
                ontology_terms = ontology_terms + convert_nones_to_list(transcript["Ontology_term"])
                dbxref_dict = parse_Dbxref(transcript["Dbxref"])
                CDS_Products = CDS_Products + transcript["CDS_Product"]
            
            # Remove dublicates among different transcripts
            ontology_terms = remove_dublicates(ontology_terms)
            CDS_Products = remove_dublicates(CDS_Products)
            dbxref_dict = remove_duplicate_values(dbxref_dict)
            
            # Write all new data (don't need the DE data) to a row
            tsv_writer.writerow([data.Genes]+
                                concat_list(gene_name)+
                                concat_list(dbxref_dict["GeneID"])+
                                concat_list(CDS_Products)+
                                concat_list(ontology_terms)+
                                concat_list(transcript_IDs)+
                                concat_list(dbxref_dict["UniProtKB"])+
                                concat_list(dbxref_dict["KEGG"])+
                                concat_list(dbxref_dict["PFAM"])+
                                concat_list(dbxref_dict["InterPro"])+
                                concat_list(dbxref_dict["EMBL"]))

if __name__ == '__main__':
    main()