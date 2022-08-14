#!/usr/bin/env python
"""
Annotate significantly differentially expressed genes using 
annotations originating from a gff3 file.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List
import pandas as pd
import numpy as np
import gffutils
from collections import namedtuple
from dataclasses import dataclass, field, asdict
from collections import Counter

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate DE genes.",
        epilog="Example: python add_gff_annotations.py ",
    )
    parser.add_argument(
        "sqlite_db",
        metavar="SQLITE_DB",
        type=Path,
        help="Sqlite DB produced with gffutils from gff3 annotations file.",
    )
    parser.add_argument(
        "de_genes",
        metavar="DE_GENES",
        type=Path,
        help="Adjusted p-value filtered DE genes data produced by deseq2.",
    )
    parser.add_argument(
        "annotated_de_genes",
        metavar="ANNOTATED_DE_GENES",
        type=Path,
        help="Annotated DE genes data in tsv file format",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def read_gff_db(db_file):
    """Read annotation sqlite database and return FeatureDB instance"""
    logger.info("Reading %s database...", db_file)

    return gffutils.FeatureDB(db_file)


def get_db_name(db_name_id: str, separator: str = ":") -> tuple:
    """
    Extract from e.g. 'FlyBase_Annotation_IDs:CG13745-PA'
    (FlyBase_Annotation_IDs,CG13745-PA)
    """
    return (db_name_id.split(separator)[0], db_name_id.split(separator)[1])


@dataclass
class Transcript:
    """Class for storing all data related to each transcript in the annotation"""

    all_dbxs = set()
    id: list[str] = field(default_factory=list)
    parent: list[str] = field(default_factory=list)
    flybase: list[str] = field(default_factory=list)
    flybase_revhits: list[str] = field(default_factory=list)
    flybase_name: list[str] = field(default_factory=list)
    gbue11: list[str] = field(default_factory=list)
    gbue11_revhits: list[str] = field(default_factory=list)
    orthodb: list[str] = field(default_factory=list)
    orthodb_revhits: list[str] = field(default_factory=list)
    ontology_term: list[str] = field(default_factory=list)
    dbxref_lst: list[str] = field(repr=False, default_factory=list)
    dbxref: dict = None

    def __post_init__(self):
        """Populate dbxref with keys and values from dbxref field in gff annotations DB"""
        if self.dbxref_lst:
            db_names: list[str] = [get_db_name(db_ref)[0] for db_ref in self.dbxref_lst]
            db_name_counts: dict = dict(Counter(db_names))
            dbxref_dict = {}
            for db_name_ref in self.dbxref_lst:
                db_name, db_id = get_db_name(db_name_ref)
                num_db_names = db_name_counts.get(db_name)
                # Add to dict if only one else update counter dict and add
                # unique ending to the dict key
                if num_db_names == 1:
                    dbxref_dict.update({db_name: db_id})
                else:
                    db_names_left: int = num_db_names - 1
                    db_name_counts.update({db_name: db_names_left})
                    new_db_name: str = f"{db_name}_{num_db_names}"
                    dbxref_dict.update({new_db_name: db_id})
            self.dbxref = dbxref_dict
            Transcript.all_dbxs.update(dbxref_dict.keys())


def create_transcripts_list(
    features_db: gffutils.FeatureDB,
) -> tuple[list[Transcript], set[str]]:
    """Iterate over all mRNA features in the features_db and create a list of transcripts.

    Args:
        features_db (gffutils.FeatureDB): FeatureDB containing all features in the gff3 database.

    Returns:
        tuple[list[Transcript], set[str]]: List of all transcripts and a set of all unique external database references.
    """
    return (
        [
            Transcript(
                id=element.attributes.get("ID"),
                parent=element.attributes.get("Parent"),
                flybase=element.attributes.get("flybase"),
                flybase_revhits=element.attributes.get("flybase_revhits"),
                flybase_name=element.attributes.get("flybase_name"),
                gbue11=element.attributes.get("gbue11"),
                gbue11_revhits=element.attributes.get("gbue11_revhits"),
                orthodb=element.attributes.get("orthodb"),
                orthodb_revhits=element.attributes.get("orthodb_revhits"),
                ontology_term=element.attributes.get("Ontology_term"),
                dbxref_lst=element.attributes.get("Dbxref"),
            )
            for element in features_db.features_of_type("mRNA")
        ],
        sorted(list(Transcript.all_dbxs)),
    )


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    sqlite_db = args.sqlite_db
    if not sqlite_db.is_file():
        logger.error("The given input file %s was not found!", sqlite_db)
        sys.exit(2)
    de_genes = args.de_genes
    if not de_genes.is_file():
        logger.error("The given input file %s was not found!", de_genes)
        sys.exit(3)

    transcripts, dbx_refs_sorted = create_transcripts_list(read_gff_db(sqlite_db))

    transcript_annotation_header: list[str] = [
        "id",
        "parent",
        "flybase",
        "flybase_revhits",
        "flybase_name",
        "gbue11",
        "gbue11_revhits",
        "orthodb",
        "orthodb_revhits",
        "ontology_term",
    ]

    df_data: list[str] = []
    header: list[str] = []

    with open(de_genes, "r") as reader:
        # Handle DE data with column names so it's easier to access each element
        diff_exp_reader = csv.reader(reader, delimiter="\t")
        de_header = next(diff_exp_reader)
        header = de_header + transcript_annotation_header + dbx_refs_sorted
        De_data = namedtuple("DE_data", de_header)
        for data in map(De_data._make, diff_exp_reader):
            # Capture all data from our current DE data row
            current_de_data = list(data)
            # Find transcripts that belong to the DE gene
            for transcript in transcripts:
                if data.Genes == transcript.parent[0]:
                    # Get each piece of transcript annotation in same order to a list
                    # so that it can be printed in correct order as annotation to the
                    # rest of DE data
                    transcript_annotation = []
                    for column in transcript_annotation_header:
                        column_value = asdict(transcript).get(column)
                        # Join all non-None column values to one string
                        if column_value:
                            transcript_annotation.append(",".join(column_value))
                        else:
                            # None:s should be written as "NA" in the annotation
                            transcript_annotation.append(np.nan)
                    # Fetch all in the DB existing dbx references from the annotation
                    for dbx_ref in dbx_refs_sorted:
                        # If we have a dbx_ref fetch values from it
                        if transcript.dbxref:
                            dbx_ref_fetched = transcript.dbxref.get(dbx_ref)
                            # If the reference value is non-None add it to the annotation
                            if dbx_ref_fetched:
                                transcript_annotation.append(dbx_ref_fetched)
                            else:
                                # None:s should be written as "NA" in the annotation
                                transcript_annotation.append(np.nan)
                        else:
                            # None:s should be written as "NA" in the annotation
                            transcript_annotation.append(np.nan)
                    # Append the data row as a list
                    df_data.append(current_de_data + transcript_annotation)

        df = pd.DataFrame(df_data, columns=header)
        # Drop all columns where all values are np.nan
        df.dropna(axis=1, how="all", inplace=True)
        # Drop duplicate "parent" column
        df.drop("parent", axis=1, inplace=True)
        # Rename columns according to below dictionary
        renamings = {
            "id": "transcript_ID",
            "flybase": "FlyBase_ID",
            "flybase_revhits": "FlyBase_reverse_hits_IDs",
            "flybase_name": "FlyBase_symbol_name",
            "ontology_term": "Annotated_Gene_Ontology_terms",
            "FlyBase": "FlyBase_reference_ID1",
            "FlyBase_2": "FlyBase_reference_ID2",
            "FlyBase_Annotation_IDs": "FlyBase_Annotation_Symbol_ID1",
            "FlyBase_Annotation_IDs_2": "FlyBase_Annotation_Symbol_ID2",
            "FlyMine": "FlyMine_ID1",
            "FlyMine_2": "FlyMine_ID2",
            "GB_protein": "GB_protein_ID1",
            "GB_protein_2": "GB_protein_ID2",
            "GB_protein_3": "GB_protein_ID3",
            "GB_protein_4": "GB_protein_ID4",
            "GB_protein_5": "GB_protein_ID5",
            "REFSEQ": "NCBI_Reference_Sequence_ID1",
            "REFSEQ_2": "NCBI_Reference_Sequence_ID2",
            "UniProt/Swiss-Prot": "UniProt_Swiss-Prot_ID",
            "UniProt/TrEMBL": "UniProt_TrEMBL_ID1",
            "UniProt/TrEMBL_2": "UniProt_TrEMBL_ID2",
        }
        df.rename(columns=renamings, inplace=True)
        df.to_csv(args.annotated_de_genes, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
