#!/usr/bin/env python
"""
Parse through a tsv file and collapse rows belonging to the same gene into one row
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from collections import namedtuple
from dataclasses import dataclass, field, asdict
from collections import Counter
from typing import Union

logger = logging.getLogger()


# def parse_args(argv=None):
#     """Define and immediately parse command line arguments."""
#     parser = argparse.ArgumentParser(
#         description="Parse through a tsv file and collapse rows belonging to the same gene into one row.",
#         epilog="Example: python each_gene_in_one_row.py combined-genes-GO_Fisher.tsv each_gene_in_one_row.tsv",
#     )
#     parser.add_argument(
#         "expanded_tsv",
#         metavar="EXPANDED_TSV",
#         type=Path,
#         help="Tsv file where data of each gene can reside in several rows.",
#     )
#     parser.add_argument(
#         "collapsed_tsv",
#         metavar="COLLAPSED_TSV",
#         type=Path,
#         help="Tsv file where each gene's data exists in one row",
#     )
#     parser.add_argument(
#         "-l",
#         "--log-level",
#         help="The desired log level (default WARNING).",
#         choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
#         default="WARNING",
#     )
#     return parser.parse_args(argv)


def open_file(filename: Path) -> list[str]:
    """Open a file into memory and return it as a list

    Args:
        filename (Path): Path to a file to be opened

    Returns:
        list[str]: List of unparsed line strings
    """
    with open(filename, "r") as f:
        return f.readlines()


def split_line(line: str, delimiter: str = "\t") -> list[str]:
    """Strip trailing new line and split string by delimiter

    Args:
        line (str): Raw line string
        delimiter (str, optional): Delimiter to split the line with. Defaults to "\t".

    Returns:
        list[str]: List of row elements
    """
    return line.strip().split(delimiter)


def parse_lines(lines: list[str], delimiter: str = "\t") -> list[list[str]]:
    """Split each unparsed string into a list and return them in a list

    Args:
        lines (list[str]): List of unparsed line strings
        delimiter (str, optional): Delimiter to split the unparsed line strings. Defaults to "\t".

    Returns:
        list[list[str]]: List of list of string row elements
    """
    return [split_line(line, delimiter) for line in lines]


def gather_unique_gene_names(
    lines: list[str], gene_name_col_index: int = 0
) -> set[str]:
    """Get unique gene names from a list of list of string row elements

    Args:
        lines (list[str]): List of list of string row elements
        gene_name_col_index (int, optional): Index position of where gene names are in the list. Defaults to 0.

    Returns:
        set[str]: Set of unique gene names
    """
    unique_gene_names: set = set()
    unique_gene_names.update([row[gene_name_col_index] for row in lines])
    return unique_gene_names


def find_num_unique(elements: tuple) -> int:
    return len(set(elements))

def remove_duplicate_elements(elements: Union[tuple,str]) -> list[str]:
    if isinstance(elements, str):
        elements: tuple = (elements,)
    return list(dict.fromkeys(elements))




def main(argv=None):
    """Coordinate argument parsing and program execution."""
    # args = parse_args(argv)
    # logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    logging.basicConfig(level="WARNING", format="[%(levelname)s] %(message)s")
    # expanded_tsv = args.expanded_tsv
    expanded_tsv = Path("combined-genes-GO_Fisher.tsv")
    if not expanded_tsv.is_file():
        logger.error("The given input file %s was not found!", expanded_tsv)
        sys.exit(2)
    lines = parse_lines(open_file(expanded_tsv))
    # Capture row header
    header = lines.pop(0)
    gene_names: set = gather_unique_gene_names(lines)
    # print(genes_data)
    # Create a dict with gene names as keys and list of row lists as values
    genes_data: dict = dict.fromkeys(gene_names, None)
    for line in lines:
        gene_name: str = line[0]
        if genes_data[gene_name]:
            genes_data[gene_name].append(line)
        else:
            genes_data[gene_name] = [line]
    
    # Collapse all data related to a gene into a dict of dicts
    for gene_data in genes_data.items():
        current_gene_name: str = gene_data[0]
        gene_data_rows: list[list[str]] = gene_data[1]
        if len(gene_data_rows) > 1:
            # Bind the column name with the 
            data_rows_packed: list[tuple[str]] = list(zip(*gene_data_rows))
            zipped: dict[dict] = dict(list(zip(header,data_rows_packed)))
            genes_data[current_gene_name] = zipped
        else:
            genes_data[current_gene_name] = dict(list(zip(header,gene_data_rows[0])))
    
    print_fields: list[str] = ["Genes","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","transcript_ID","FlyBase_ID","FlyBase_reverse_hits_IDs","FlyBase_symbol_name","gbue11","gbue11_revhits","orthodb","orthodb_revhits","Annotated_Gene_Ontology_terms","FlyBase_reference_ID1","FlyBase_reference_ID2","FlyBase_Annotation_Symbol_ID1","FlyBase_Annotation_Symbol_ID2","FlyMine_ID1","FlyMine_ID2","GB_protein_ID1","GB_protein_ID2","GB_protein_ID3","GB_protein_ID4","GB_protein_ID5","NCBI_Reference_Sequence_ID1","NCBI_Reference_Sequence_ID2","UniProt_Swiss-Prot_ID","UniProt_TrEMBL_ID1","UniProt_TrEMBL_ID2","modMine","modMine_2","GO_ID_BP","Term_BP","Annotated_BP","Significant_BP","Expected_BP","Rank_in_classicFisher_BP","classicFisher_BP","elimFisher_BP","weight01Fisher_BP","parentchildFisher_BP","weightFisher_BP","leaFisher_BP","GO_ID_CC","Term_CC","Annotated_CC","Significant_CC","Expected_CC","Rank_in_classicFisher_CC","classicFisher_CC","elimFisher_CC","weight01Fisher_CC","parentchildFisher_CC","weightFisher_CC","leaFisher_CC","GO_ID_MF","Term_MF","Annotated_MF","Significant_MF","Expected_MF","Rank_in_classicFisher_MF","classicFisher_MF","elimFisher_MF","weight01Fisher_MF","parentchildFisher_MF","weightFisher_MF","leaFisher_MF"]
    
    for data in genes_data.items():
        current_gene_name: str = data[0]
        gene_data: dict = data[1]
        deduped = []
        for data_field in gene_data.items():
            field_header: str = data_field[0]
            data_objects: tuple = data_field[1]
            # if field_header in shortenable_fields:
            deduped: list[str] = remove_duplicate_elements(data_objects)
            # else:
            #     deduped: list[str] = list(set(data_objects))
            genes_data[current_gene_name][field_header] = deduped
    
    with open('test.tsv', 'w') as f:
        # create the csv writer
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(print_fields)
        for data in genes_data.items():
            current_gene_name: str = data[0]
            gene_data: dict = data[1]
            # [";".join(gene_data.get(data_value)) for data_value in print_fields]
            # write a row to the csv file
            writer.writerow([";".join(gene_data.get(data_value)) for data_value in print_fields])
    
    print("test")
if __name__ == "__main__":
    sys.exit(main())
