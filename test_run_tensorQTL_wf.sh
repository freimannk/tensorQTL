#!/bin/bash

#Töö nimi
#SBATCH -J GTEx

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 11:40:00

# mälu
#SBATCH --mem=5G

module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load squashfs/4.4


NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/LCL/GTEx/GTEx_naive_LCL --vcf /gpfs/space/home/a82371/LCL/merged_freq_filtered_extracted_EUR_bial --covariates /gpfs/space/home/a82371/LCL/GTEx/GTEx_LCL_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/qcnorm/LCL_naive_results/GTEx/GTEx_LCL/normalised/ge/qtl_group_split_norm/CAP.LCL_naive.tsv --study GTEx_naive_LCL --sample_genotype_ids /gpfs/space/home/a82371/LCL/GTEx/LCL_GTEx_study_ids.txt --only_autosomal_chr false --median_tpm_filtration_file /gpfs/space/home/a82371/LCL_MULTI/LCL_naive_median_tpm_filtered_genes.tsv.gz -profile tartu_hpc -resume
