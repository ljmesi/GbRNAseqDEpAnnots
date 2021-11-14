#!/usr/bin/env python3
"""
Read gff-file and output a qffutils usable sqlite database
"""

import click
import gffutils

def create_db(gff_file, db_file):
    """Read a gff-file and create FeatureDB sqlite db instance"""
    click.echo(f"Creating {db_file} database")
    gffutils.create_db(gff_file, 
                        dbfn=db_file, 
                        force=True, 
                        keep_order=True,
                        merge_strategy='merge', 
                        sort_attribute_values=True)
    click.echo(f"{db_file} database created")

@click.command()
@click.argument('gff-file', type=click.Path(exists=True))
@click.argument('output-db-file', type=click.Path())
def main(gff_file, output_db_file):
    """Run the command line program"""
    create_db(gff_file, output_db_file)

if __name__ == '__main__':
    main()