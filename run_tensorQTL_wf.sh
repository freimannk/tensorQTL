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

NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/BLUEPRINT/Blueprint_naive_Tcells_144509572_v2 --vcf /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT.MAF001 --covariates /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_T_cells.covariates_added_sex.txt --expression_file /gpfs/space/home/a82371/qcnorm/results_rnaseq_naive_Tcells/Blueprint/Blueprint_naive_Tcells/normalised/ge/qtl_group_split_norm/BLUEPRINT_PE.T-cell.tsv --study Blueprint_T_cell --sample_genotype_ids /gpfs/space/home/a82371/BLUEPRINT/BLUEPRINT_naive_Tcells_ids.txt --only_autosomal_chr true --median_tpm_filtration_file /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/T-cell_median_tpm_filtered_genes.tsv.gz -profile tartu_hpc -resume
