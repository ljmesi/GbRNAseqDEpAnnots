#!/usr/bin/env python3
"""
Read gff-file and output a qffutils usable sqlite database
"""

import click

def create_GO_set(genes_GOs_file, GO_term_column):
    """Collect all unique GO-terms in a set"""
    GO_set = set()
    for line in genes_GOs_file:
        line_list = line.strip().split("\t")
        # Check if there aren't any GO-terms for this gene
        if len(line_list) > GO_term_column and line_list[GO_term_column]:
            # Parse comma separated GO-terms
            GO_list = line_list[GO_term_column].split(",")
            for GO_term in GO_list:
                #print(GO_term)
                GO_set.add(GO_term)
    return GO_set

def print_GO_set(GO_set, outfile):
    """Iterate over all unique GO-terms and print them to out-file"""
    for GO_term in GO_set:
        outfile.write(GO_term + "\n")

@click.command()
@click.argument('genes-gos-file', type=click.File('r'))
@click.argument('output-file', type=click.File('w'))
@click.option('-c',
            '--go-term-column',
            type=int,
            default=1,
            help="Column where comma separated GO-terms are (should be the last one)",
            required=True)
def main(genes_gos_file, output_file, go_term_column):
    """Run the command line program"""
    GO_set = create_GO_set(genes_gos_file, go_term_column)
    print_GO_set(GO_set, output_file)

if __name__ == '__main__':
    main()
