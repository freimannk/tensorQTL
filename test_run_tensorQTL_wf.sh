#!/bin/bash

#Töö nimi
#SBATCH -J GTEx

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 11:40:00

# mälu
#SBATCH --mem=5G

#siia alla on vaja kirjutada protsessid
module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load squashfs/4.4

#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath  --vcf  --covariates  --expression_file  --study  --sample_genotype_ids  --only_autosomal_chr true -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/BLUEPRINT/naive_Tcells --vcf /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT.MAF001 --covariates /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE_ge_T-cell.covariates.txt --expression_file /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE.T-cell.tsv --study BLUEPRINT_naive_Tcells --sample_genotype_ids /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_naive_Tcells_ids.txt --only_autosomal_chr true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz --filter_cis_window 5000000 -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/BLUEPRINT/naive_Tcells_1000var --vcf /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT.MAF001 --covariates /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE_ge_T-cell.covariates.txt --expression_file /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE.T-cell.tsv --study BLUEPRINT_naive_Tcells --sample_genotype_ids /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_naive_Tcells_ids.txt --only_autosomal_chr true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz --filter_cis_window 5000000 -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/CEDAR/naive_CD4_Tcells --vcf /gpfs/space/home/a82371/CEDAR/CEDAR.filtered --covariates /gpfs/space/home/a82371/CEDAR/CEDAR_microarray_T-cell_CD4.covariates.txt --expression_file /gpfs/space/home/a82371/CEDAR/CEDAR_T-cell_CD4_gene_expression.tsv --study CEDAR_T-cell_CD4 --sample_genotype_ids /gpfs/space/home/a82371/CEDAR/CEDAR_T-cell_CD4_ids.txt --only_autosomal_chr false --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz --filter_cis_window 5000000 -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/INTERVAL/naive_plasma --vcf /gpfs/space/home/a82371/INTERVAL/INTERVAL.filtered --covariates /gpfs/space/home/a82371/INTERVAL/INTERVAL.covariates.txt --expression_file /gpfs/space/home/a82371/INTERVAL/INTERVAL_protein_expression.tsv --study INTERVAL_naive_plasma --sample_genotype_ids /gpfs/space/home/a82371/INTERVAL/INTERVAL_plasma_protein_ids.txt --only_autosomal_chr true --filter_cis_window 5000000 -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/LCL_MULTI/naive_LCL_test_cis_filt2 --vcf /gpfs/space/home/a82371/LCL_MULTI/merged_freq_filtered_extracted_EUR_bial --covariates /gpfs/space/home/a82371/LCL_MULTI/LCL_multi_ge_naive.covariates.txt --expression_file /gpfs/space/home/a82371/LCL_MULTI/LCL_multi_ge_naive.tsv --study naive_LCL_MULTI --sample_genotype_ids /gpfs/space/home/a82371/LCL_MULTI/LCL_multy_study_ids.txt --only_autosomal_chr true --missing_DS true --median_tpm_filtration_file /gpfs/space/home/a82371/LCL_MULTI/LCL_naive_median_tpm_filtered_genes.tsv.gz --filter_cis_window 5000000 -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/BLUEPRINT/Blueprint_naive_Tcells --vcf /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT.MAF001 --covariates /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_T_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/qcnorm/results_rnaseq_naive_Tcells/Blueprint/Blueprint_naive_Tcells/normalised/ge/qtl_group_split_norm/BLUEPRINT_PE.T-cell.tsv --study Blueprint_T_cell --sample_genotype_ids /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_naive_Tcells_ids.txt --only_autosomal_chr true --missing_DS true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/CEDAR/CEDAR_naive_CD4_Tcells --vcf /gpfs/space/home/a82371/CEDAR/CEDAR.filtered --covariates /gpfs/space/home/a82371/CEDAR/CEDAR_cd4_T_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/CEDAR/CEDAR.T-cell_CD4.HumanHT-12_V4_norm_exprs_translated.tsv --study CEDAR_T-cell_CD4 --sample_genotype_ids /gpfs/space/home/a82371/CEDAR/CEDAR_T-cell_CD4_ids.txt --only_autosomal_chr false --missing_DS true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz -profile tartu_hpc -resume
#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/LCL/GTEx/GTEx_naive_LCL --vcf /gpfs/space/home/a82371/LCL/merged_freq_filtered_extracted_EUR_bial --covariates /gpfs/space/home/a82371/LCL/GTEx/GTEx_LCL_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/qcnorm/LCL_naive_results/GTEx/GTEx_LCL/normalised/ge/qtl_group_split_norm/CAP.LCL_naive.tsv --study GTEx_naive_LCL --sample_genotype_ids /gpfs/space/home/a82371/LCL/GTEx/LCL_GTEx_study_ids.txt --only_autosomal_chr false --missing_DS true --median_tpm_filtration_file /gpfs/space/home/a82371/LCL_MULTI/LCL_naive_median_tpm_filtered_genes.tsv.gz --maf_filter 0 -profile tartu_hpc -resume
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/LCL/GTEx/GTEx_naive_LCL --vcf /gpfs/space/home/a82371/LCL/merged_freq_filtered_extracted_EUR_bial --covariates /gpfs/space/home/a82371/LCL/GTEx/GTEx_LCL_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/qcnorm/LCL_naive_results/GTEx/GTEx_LCL/normalised/ge/qtl_group_split_norm/CAP.LCL_naive.tsv --study GTEx_naive_LCL --sample_genotype_ids /gpfs/space/home/a82371/LCL/GTEx/LCL_GTEx_study_ids.txt --only_autosomal_chr false --missing_DS true --median_tpm_filtration_file /gpfs/space/home/a82371/LCL_MULTI/LCL_naive_median_tpm_filtered_genes.tsv.gz -profile tartu_hpc -resume
