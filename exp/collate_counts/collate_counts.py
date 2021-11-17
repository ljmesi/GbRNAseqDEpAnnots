#!/usr/bin/env python3
"""
Read raw STAR-aligner gene counts tsv files and 
concatenate gene counts to one output tsv file.
"""

import argparse
import pandas as pd
from pathlib import Path

def write_tsv(filename, df):
    """Write pandas DataFrame as tsv file"""
    df.to_csv(filename, sep='\t')

def get_file_stem(file_name):
    """Extract as string file stem name from file name"""
    return str(Path(file_name).stem)

def parse_sample_name(file_name):
    """Obtain sample name from file name"""
    file_stem = get_file_stem(file_name)    
    slist = file_stem.split('_')[0].split('-')
    return slist[1] + "-" + slist[2]

def read_tsv(filename, header_rows_count):
    """Read tsv file and return pandas DataFrame"""
    sample_name = parse_sample_name(filename)
    df = pd.read_csv(filename, 
                    sep = '\t',
                    usecols = [0,1],
                    names = ["-",sample_name],
                    skiprows = header_rows_count,
                    index_col = 0)
    df.index.names = ['Genes']
    return df

def concat_cols(df_list):
    """Concatenate list of one column pandas DataFrames with same index"""
    return pd.concat(df_list, axis=1)

def order_cols(df, col_name_list):
    """Rearrange pandas DataFrame columns in the same order as given column name list"""
    return df[col_name_list]

def main(args):
    """Run the command line program"""
    unordered_df = concat_cols([read_tsv(file_path, 4) for file_path in args.infiles])
    df = order_cols(unordered_df, ["12-5","12-6","12-13","12-15","18-1","18-2","18-3","18-8"])
    write_tsv(args.output, df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infiles',
                        nargs='*', 
                        help='Input file names')
    parser.add_argument('-o', '--output',
                        default=None,
                        help='Output tsv file name')
    args = parser.parse_args()
    main(args)
