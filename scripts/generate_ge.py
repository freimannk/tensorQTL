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


def put_chr_in_front(chr_nr):
    return 'chr' + str(chr_nr)


chr = merged_df.chr.map(put_chr_in_front)
tensorQTL_ge_df = merged_df.assign(chr=chr.values)
tensorQTL_ge_df = tensorQTL_ge_df.rename(columns={'chr': '#chr'}, inplace=False)
tensorQTL_ge_df.to_csv(f"{args.study_name}_tensorQTL_wf_ge.bed", sep="\t", quoting=csv.QUOTE_NONE,
                       index=False)
