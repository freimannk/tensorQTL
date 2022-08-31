#!/usr/bin/env python3

import argparse
import csv
import pandas as pd

parser = argparse.ArgumentParser(description="Generate gene expression file with correct format.")

parser.add_argument('-g', '--ge_file', required=True, type=str,
                    help="gene expression input file.")
parser.add_argument('-s', '--samples_genotype_ids', required=True, type=str,
                    help="Sample genotype_ids file.")
parser.add_argument('-t', '--genes_tss_file', required=True, type=str,
                    help="Homo sapiens GRCh38 96 genes TSS")
parser.add_argument('-n', '--study_name', required=True, type=str,
                    help="Study name")

parser.add_argument('-f', '--median_tpm_filtration_file', nargs='?', const=None, type=str,
                    help="median_tpm_filtration_file path to filter genes based on tpm")

args = parser.parse_args()

samples_genotype_ids_df = pd.read_csv(args.samples_genotype_ids, sep='\t', index_col=0)

ge_df = pd.read_csv(args.ge_file, sep='\t', index_col=0)
samples_to_include = samples_genotype_ids_df.index.tolist()
ge_df = ge_df[samples_to_include]
samples_id_genotype_id_dic = pd.Series(samples_genotype_ids_df.genotype_id.values,
                                       index=samples_genotype_ids_df.index).to_dict()
ge_df = ge_df.rename(columns=samples_id_genotype_id_dic, inplace=False)

Homo_sapiens_GRCh38_96_genes_TSS_df = pd.read_csv(args.genes_tss_file, sep='\t', index_col=0)

merged_df = pd.merge(Homo_sapiens_GRCh38_96_genes_TSS_df,
                     ge_df,
                     right_index=True,
                     how='inner',
                     left_index=True
                     )

tensorQTL_ge_df = merged_df.rename(columns={'chr': '#chr'}, inplace=False)

if args.median_tpm_filtration_file is not None:
    tpm_filtered_genes_df = pd.read_csv(args.median_tpm_filtration_file, sep='\t', compression='gzip')
    phenotype_ids = tpm_filtered_genes_df.phenotype_id.to_list()
    tensorQTL_ge_df = tensorQTL_ge_df[tensorQTL_ge_df.gene_id.isin(phenotype_ids)]

tensorQTL_ge_df= tensorQTL_ge_df.sort_values(by=['#chr', 'start'])
tensorQTL_ge_df.to_csv(f"{args.study_name}_tensorQTL_wf_ge.bed", sep="\t", quoting=csv.QUOTE_NONE,
                       index=False)