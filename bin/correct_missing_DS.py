#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import argparse

parser = argparse.ArgumentParser(description="Input file list and variants pattern.")

parser.add_argument('-i', '--inputFile', required=True,
                    help="Input file path.")
parser.add_argument('-n', '--fileName', required=True, type=str,
                    help="File name.")

args = parser.parse_args()



df = pd.read_csv(args.inputFile, sep='\t', compression='gzip',low_memory=False)
df_w_nans = df.applymap(lambda x: np.nan if x == '.' else x)
df_w_nans.set_index(['CHROM','POS', 'ID'], inplace=True)
df_w_avg_DS = df_w_nans.T.fillna(df_w_nans.mean(axis=1)).T
df_w_avg_DS = df_w_avg_DS.reset_index()
df_w_avg_DS.to_csv(f"{args.fileName}.txt.gz", sep="\t", quoting=csv.QUOTE_NONE,compression='gzip',
                         index=False)