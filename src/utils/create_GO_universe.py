#!/usr/bin/env python3
"""
Read gff-file and output a GO term universe as tsv file
"""

import click
import pandas as pd
from pathlib import Path
import gffutils

def read_gff_db(db_file):
    """Read annotation sqlite database and return FeatureDB instance"""
    click.echo(f"Reading {db_file} database...")
    return gffutils.FeatureDB(db_file)

def remove_dublicates(list_of_GOs):
    """Remove dublicate GO terms from a list"""
    #click.echo("Removing duplicate GO-terms...")
    return list(set(list_of_GOs))

def join_dict_value_lists(all_genes):
    """
    Iterate over all list of lists in dict values and update the value to
    one concatenated list
    """
    click.echo("Joining list of duplicated GO-terms into one and removing duplicates...")
    for key,value in all_genes.items():
        if value is not None:
            all_genes[key] = remove_dublicates(sum(value, []))
    click.echo("Lists concatenated and GO-terms deduplicated")
    return all_genes

def collect_all_mRNAs(db, gene):
    """Gather all mRNAs of a gene into lists and return them as a list"""
    #click.echo("Collecting all mRNAs from genes...")
    mRNAs = []
    for feature in db.children(gene):
        if(feature.featuretype == "mRNA"):
            mRNAs.append(feature)
    #click.echo("All mRNAs from genes collected")
    return mRNAs

def populate_GO_dict(db):
    """
    Read FeatureDB and return a dict of gene IDs as keys and unique GO terms as
    values
    """
    click.echo("Populating GO-terms dictionary...")
    all_genes = {}
    for gene in db.features_of_type('gene'):
        for mRNA in collect_all_mRNAs(db, gene):
            GO_terms = mRNA.attributes.get("Ontology_term")
            if GO_terms is None:
                GO_terms = []
            if gene.id in all_genes:
                all_genes[gene.id].append(GO_terms)
            else:
                all_genes[gene.id] = [GO_terms]
    click.echo("GO-terms dictionary populated")
    return join_dict_value_lists(all_genes)
    
# def fill_Nones(GO_dict, fill_value):
#     """"""
#     for key in GO_dict:
#         if not GO_dict[key]:
#             GO_dict[key] = fill_value
#     return GO_dict

def concatenate_unique_GO_list(GO_dict):
    """
    Concatenate all GO-terms in a list to one string 
    and return the updated dict
    """
    click.echo("Joining together GO-terms in a list into one string...")
    for key in GO_dict:
        if isinstance(GO_dict[key], list):
            GO_dict[key] = ",".join(GO_dict[key])
    click.echo("GO-terms list concatenated")
    return GO_dict

def convert_dict_to_df(GO_dict):
    """Convert each dict entry to one row in pandas df and return the whole df"""
    click.echo("Converting GO-terms dict to pandas dataframe...")
    df = pd.DataFrame.from_dict(GO_dict, 
                                    orient='index',
                                    columns=['GO-terms'])
    df.index.name = "transcripts"
    click.echo("GO-terms dict to pandas dataframe converted")
    return df

def write_tsv(filename, df):
    """Write pandas DataFrame as tsv file"""
    click.echo(f"Writing into {filename} file pandas dataframe...")
    df.to_csv(filename, sep='\t')
    click.echo(f"Pandas dataframe writen to {filename}")

@click.command()
@click.argument('db-file', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(db_file, output):
    """Run the command line program"""
    db = read_gff_db(db_file)
    click.echo(f"{db_file} database read")
    GO_dict = populate_GO_dict(db)
    click.echo("GO-terms dictionary populated, joined and deduped")
    GO_dict = concatenate_unique_GO_list(GO_dict)
    GO_df = convert_dict_to_df(GO_dict)
    write_tsv(output, GO_df)

if __name__ == '__main__':
    main()