#!/usr/bin/env python3
"""
Read DE results tsv file and add annotations from annotations sqlite database
"""

import click
from pathlib import Path
import csv
import sys
from collections import namedtuple

def open_file(filename):
    """Open a file into memory and return it as a list"""
    with open(filename, 'r') as f:
        return f.readlines()

def parse_line(line, delimiter = '\t'):
    return line.strip().split(delimiter)

def parse_lines(lines, delimiter = '\t'):
    """"""
    return [parse_line(line, delimiter) for line in lines]

def pop_parsed_header(lines, delimiter = '\t'):
    return (lines, parse_line(lines.pop(0), delimiter))

def fetch_nth_elements(lines_ll, n):
    """Return a list of n:th elements in a lines list of lists"""
    return [line[n] for line in lines_ll]

def add_to_genes_dict(genes, gene_name, gene_data):
    """Add to a genes dictionary genes data lists"""
    if gene_name in genes:
        genes[gene_name].append(gene_data)
    else:
        genes[gene_name] = [gene_data]
    return genes

def build_genes_dict(lines_ll):
    """Create a genes dictionary with a genes data list of lists"""
    genes = {}
    for line in lines_ll:
        add_to_genes_dict(genes, line.pop(0), line)
    return genes

def create_gene_info_lists(info_ll):
    # Create a list of dictionaries for keeping unique row values in the same order
    info_list = []
    # Create as many empty dictionaries as there are elements in the first list
    for index in info_ll[0]:
        info_list.append(dict())
    
    # Populate each dictionary with same row position data
    for info_l in info_ll:
        for idx, data in enumerate(info_l):
            info_list[idx].update({data: None})
    return info_list

def create_ordered_gene_infos(genes):
    gene_infos = {}
    for gene, info_ll in genes.items():
        gene_infos.update({gene: create_gene_info_lists(info_ll)})
    return gene_infos

def concatenate_dict_keys(data_dict, delimiter = ','):
    return delimiter.join(data_dict.keys())

def convert_dicts_to_strings(info_l):
    return [concatenate_dict_keys(info) for info in info_l]

def simplify_nonredundant_gene_infos(genes):
    nonredundant_gene_infos = {}
    for gene, info_l in genes.items():
        nonredundant_gene_infos.update({gene: convert_dicts_to_strings(info_l)})
    return nonredundant_gene_infos

def write_tsv(filename, original_header, genes_data, out_file_header = ['Genes', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'Gene_name', 'GeneID', 'CDS_Products', 'Gff_ontology_terms', 'Transcript_IDs', 'UniProtKB_Swiss_Prot', 'KEGG', 'PFAM', 'InterPro', 'EMBL', 'GO_ID_BP', 'Term_BP', 'GO_ID_CC', 'Term_CC', 'GO_ID_MF', 'Term_MF']):
    # Use namedtuples so output file columns can be chosen by their names
    Data = namedtuple("Data", original_header)
    with open(filename, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        # Write header
        tsv_writer.writerow(out_file_header)
        for gene_name, data in genes_data.items():
            data_nt = Data._make([gene_name]+data)
            tsv_writer.writerow([gene_name,
                                data_nt.baseMean, 
                                data_nt.log2FoldChange, 
                                data_nt.lfcSE, 
                                data_nt.stat, 
                                data_nt.pvalue, 
                                data_nt.padj, 
                                data_nt.Gene_name, 
                                data_nt.GeneID, 
                                data_nt.CDS_Products, 
                                data_nt.Gff_ontology_terms, 
                                data_nt.Transcript_IDs, 
                                data_nt.UniProtKB_Swiss_Prot, 
                                data_nt.KEGG, 
                                data_nt.PFAM, 
                                data_nt.InterPro, 
                                data_nt.EMBL, 
                                data_nt.GO_ID_BP, 
                                data_nt.Term_BP, 
                                data_nt.GO_ID_CC, 
                                data_nt.Term_CC, 
                                data_nt.GO_ID_MF, 
                                data_nt.Term_MF])

@click.command()
@click.argument("in-file", type=click.Path(exists=True))
@click.argument('output-file', type=click.Path())
def main(in_file, output_file):
    """Run the main CLI application"""
    lines = open_file(in_file)
    lines, header = pop_parsed_header(lines)
    lines = parse_lines(lines)
    genes = build_genes_dict(lines)
    nonredundant_genes_dict = create_ordered_gene_infos(genes)
    genes_data = simplify_nonredundant_gene_infos(nonredundant_genes_dict)
    write_tsv(output_file,header,genes_data)

if __name__ == '__main__':
    main()