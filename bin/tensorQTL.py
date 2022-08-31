#!/usr/bin/env python3

import argparse
import os

import pandas as pd
import tensorqtl
import torch
from tensorqtl import trans

print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

parser = argparse.ArgumentParser(description="Input file list and variants pattern.")

parser.add_argument('-e', '--expression_file', required=True,
                    help="Expression file.")
parser.add_argument('-c', '--covariates_file', required=True,
                    help="Expression file.")
parser.add_argument('-g', '--genotype_matrix_file', required=True, type=str,
                    help="Genotype matrix file (GT or GS).")
parser.add_argument('-p', '--pvalue', required=True, type=float,
                    help="Pvalue to filter.")
parser.add_argument('-n', '--study_name', required=True, type=str,
                    help="Study name.")
parser.add_argument('-f', '--filter_cis_window', nargs='?', const=None, type=int,
                    help="cis window to filter out from results.")
parser.add_argument('-m', '--maf', required=True, type=float,
                    help="MAF to filter.")

args = parser.parse_args()
df = pd.read_csv(args.genotype_matrix_file, sep='\t', compression='gzip')
#df['CHROM'] = 'chr' + df['CHROM'].astype(str)
df = df.rename(columns={'POS': 'pos', 'ID': 'snp', 'CHROM': 'chrom'})
df = df.set_index('snp')
variant_df = df[['chrom', 'pos']]
df = df.rename_axis('iid', axis=1)
genotype_df = df.drop(columns=['chrom', 'pos'])





if len(genotype_df.index) > 0:
# load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.expression_file)
    covariates_df = pd.read_csv(args.covariates_file, sep='\t', index_col=0).T
    trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
                            return_sparse=True, pval_threshold=args.pvalue, maf_threshold=args.maf)
    if args.filter_cis_window is not None:
        trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=args.filter_cis_window)

    trans_df.to_parquet(os.path.join('.', args.study_name + '.trans_qtl_pairs.parquet'))

