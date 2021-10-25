#!/bin/bash

#Töö nimi
#SBATCH -J BLUEPRINT_Tcells

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 40:40:00

# mälu
#SBATCH --mem=1G

#siia alla on vaja kirjutada protsessid


module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4




#NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath  --vcf  --covariates  --expression_file  --study  --sample_genotype_ids  --only_autosomal_chr true -profile tartu_hpc -resume
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/BLUEPRINT/naive_Tcells --vcf /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT.MAF001 --covariates /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE_ge_T-cell.covariates.txt --expression_file /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_PE.T-cell.tsv --study BLUEPRINT_naive_Tcells --sample_genotype_ids /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_naive_Tcells_ids.txt --only_autosomal_chr true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz --filter_cis_window 5000000 -profile tartu_hpc -resume
