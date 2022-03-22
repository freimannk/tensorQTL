import pandas as pd
import argparse
import csv


parser = argparse.ArgumentParser(description="Generate gene expression file with correct format.")


parser.add_argument('-s', '--samples_genotype_ids', required=True, type=str,
                    help="Sample genotype_ids file.")
parser.add_argument('-n', '--study_name', required=True, type=str,
                    help="Study name.")
parser.add_argument('-c', '--covariate_file', nargs='?', const=None, type=str,
                    help="Covariates file.")

args = parser.parse_args()
cov_df = pd.read_csv(args.covariate_file, sep='\t', index_col=0)
samples_genotype_ids_df = pd.read_csv(args.samples_genotype_ids, sep='\t', index_col=0)
samples_to_include = samples_genotype_ids_df.genotype_id.tolist()
cov_filtered_df = cov_df[samples_to_include]
cov_filtered_df.to_csv(f"{args.study_name}_filtered_covariates.txt", sep="\t", quoting=csv.QUOTE_NONE,
                       index=True)