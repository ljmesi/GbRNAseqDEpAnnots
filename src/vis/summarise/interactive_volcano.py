#!/usr/bin/env python
"""
Create an interactive volcano plot from differential expression data and
annotatations. 
"""


import pandas as pd
import math
import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import plotly.graph_objects as go

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create an interactive volcano plot.",
        epilog="Example: python interactive_volcano.py annot.tsv de_data.tsv volcano.html",
    )
    parser.add_argument(
        "annotations",
        metavar="ANNOTATIONS",
        type=Path,
        help="Annotations tsv file.",
    )
    parser.add_argument(
        "de_data",
        metavar="DE_DATA",
        type=Path,
        help="All (i.e. not filtered by adjusted p-value) DE genes data produced by deseq2.",
    )
    parser.add_argument(
        "volcano_plot",
        metavar="VOLCANO_PLOT",
        type=Path,
        help="Volcano plot which includes annotations for de genes",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    annot = args.annotations
    if not annot.is_file():
        logger.error("The given input file %s was not found!", annot)
        sys.exit(2)
    de_data = args.de_data
    if not de_data.is_file():
        logger.error("The given input file %s was not found!", de_data)
        sys.exit(3)

    df = pd.read_table(annot, index_col=False, dtype={"padj": float})
    df["Significance"] = "-log_10(padj)<=0.01"

    df_all = pd.read_table(de_data, index_col=False, dtype={"padj": float})
    df_all["Significance"] = "Not significant"

    data = df_all.merge(df, how="left", on="Genes", suffixes=["", "_y"], copy=True)

    def significator(row_val):
        def isNaN(num):
            return num != num

        if not isNaN(row_val):
            return "-log_10(padj)<=0.01"
        return "Not significant"

    data["Significance"] = data.apply(
        lambda row: significator(row.Significance_y), axis=1
    )

    # Drop duplicate columns
    cols = [
        "baseMean_y",
        "log2FoldChange_y",
        "lfcSE_y",
        "stat_y",
        "pvalue_y",
        "Significance_y",
        "padj_y",
    ]
    data.drop(cols, inplace=True, axis=1)

    data["-log_10_padj"] = data.apply(lambda row: math.log10(row.padj) * -1, axis=1)

    nanable_cols = [
        "baseMean",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ]

    for col in nanable_cols:
        data.loc[data["Significance"] == "Not significant", col] = np.nan
    data.fillna("", inplace=True)

    hover_cols = {
        "Genes": "",
        "baseMean": "Base mean",
        "log2FoldChange": "Log2 fold change",
        "lfcSE": "Log fold change SE",
        "stat": "Stat",
        "pvalue": "P-value",
        "padj": "Adj p-value",
        "-log_10_padj": "-log10 adj p-value",
        "transcript_ID": "Transcript ID",
        "FlyBase_ID": "FlyBase ID",
        "FlyBase_reverse_hits_IDs": "FlyBase reverse hits IDs",
        "FlyBase_symbol_name": "FlyBase symbol name",
        "gbue11": "gbue11",
        "gbue11_revhits": "gbue11 revhits",
        "orthodb": "Orthodb",
        "orthodb_revhits": "Orthodb revhits",
        "Annotated_Gene_Ontology_terms": "Annotated Gene Ontology terms",
        "FlyBase_reference_ID1": "FlyBase reference ID 1",
        "FlyBase_reference_ID2": "FlyBase reference ID 2",
        "FlyBase_Annotation_Symbol_ID1": "FlyBase Annotation Symbol ID 1",
        "FlyBase_Annotation_Symbol_ID2": "FlyBase Annotation Symbol ID 2",
        "FlyMine_ID1": "FlyMine ID 1",
        "FlyMine_ID2": "FlyMine ID 2",
        "GB_protein_ID1": "GB protein ID 1",
        "GB_protein_ID2": "GB protein ID 2",
        "GB_protein_ID3": "GB protein ID 3",
        "GB_protein_ID4": "GB protein ID 4",
        "GB_protein_ID5": "GB protein ID 5",
        "NCBI_Reference_Sequence_ID1": "NCBI Reference Sequence ID 1",
        "NCBI_Reference_Sequence_ID2": "NCBI Reference Sequence ID 2",
        "UniProt_Swiss-Prot_ID": "UniProt Swiss-Prot ID",
        "UniProt_TrEMBL_ID1": "UniProt TrEMBL ID 1",
        "UniProt_TrEMBL_ID2": "UniProt TrEMBL ID 2",
        "modMine": "ModMine",
        "modMine_2": "ModMine 2",
        "GO_ID_BP": "BP GO IDs",
        "Term_BP": "BP terms",
        "GO_ID_CC": "CC GO IDs",
        "Term_CC": "CC terms",
        "GO_ID_MF": "MF GO IDs",
        "Term_MF": "MF terms",
    }

    hover_texts: list[str] = []
    for index, row in data.iterrows():
        # For each row we want a fresh list to keep track of what values
        # are in each column of the row
        hover_texts_row: list[str] = []
        for index, element in enumerate(row):
            col_name: str = row.index[index]  # e.g. baseMean
            row_value: str = element  # e.g. 50282.2933422622
            hr_colname = hover_cols.get(col_name)  # e.g. Base mean
            # Present gene names as headers bold face headers in the hover text
            if col_name == "Genes":
                hover_texts_row.append(f"<b>{row_value}</b><br>")
                continue
            # Add a hover text only if the row value is among the wished hover text columns
            # and if there is a something to show in the first place
            if hr_colname and row_value:
                hover_texts_row.append(f"{hr_colname}: {row_value}<br>")
        # Produce the hover text from its elements
        hover_text = "".join(hover_texts_row)
        # Remove the last <br> from the hover text
        hover_text = hover_text[: len(hover_text) - 4]
        hover_texts.append(hover_text)

    # Add the hover texts as a column to the whole df so they go along
    data.insert(2, "hover_text", hover_texts)

    df_not_sig = data[data["Significance"] == "Not significant"]
    df_sig = data[data["Significance"] == "-log_10(padj)<=0.01"]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df_sig["log2FoldChange"],
            y=df_sig["-log_10_padj"],
            mode="markers",
            hovertext=df_sig["hover_text"].tolist(),
            hoverinfo="text",
            name="padj<=0.01",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=df_not_sig["log2FoldChange"],
            y=df_not_sig["-log_10_padj"],
            mode="markers",
            hovertext=df_not_sig["hover_text"].tolist(),
            hoverinfo="text",
            name="Not significant",
        )
    )
    fig.update_layout(
        title_text="Volcano plot",
        template="plotly_white",
        title="Volcano plot",
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text="Log2 fold change<br><sup>(Genes with positive value are upregulated in 18 hours light condition)</sup>"
            )
        ),
        yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="-log10(adj p-value)")),
    )
    fig.write_html(args.volcano_plot, config={"displaylogo": False})


if __name__ == "__main__":
    sys.exit(main())
