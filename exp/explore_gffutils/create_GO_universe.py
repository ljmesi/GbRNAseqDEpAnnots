#!/usr/bin/env python3
"""
Read gff-file and output a GO term universe as tsv file
"""

import argparse
import pandas as pd
from pathlib import Path
import gffutils


def create_db(gff_file):
    return gffutils.create_db(gff_file, 
                                dbfn='Gb.db', 
                                force=True, 
                                keep_order=True,
                                merge_strategy='merge', 
                                sort_attribute_values=True)

def populate_GO_dict(db):
    all_mRNAs = {}
    for mRNA in db.features_of_type('mRNA'):
        all_mRNAs[mRNA.id] = mRNA.attributes.get("Ontology_term")
    return all_mRNAs

def convert_dict_to_df(GO_dict):
    return pd.DataFrame.from_dict(GO_dict, 
                                    orient='index',
                                    columns=['transcript', 'GO-term'])

def write_tsv(filename, df):
    """Write pandas DataFrame as tsv file"""
    df.to_csv(filename, sep='\t')


def main(args):
    """Run the command line program"""
    write_tsv(
        args.output, 
        convert_dict_to_df(
            populate_GO_dict(
                create_db(args.gff)
                )
            )
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gff',
                        type=argparse.FileType('r'),
                        nargs='*', 
                        default='-', 
                        help='Input gff-file name')
    parser.add_argument('-o', '--output',
                        default=None,
                        help='Output GO term universe tsv file name')
    args = parser.parse_args()
    main(args)