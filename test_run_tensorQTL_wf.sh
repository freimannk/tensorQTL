#!/bin/bash

#Töö nimi
#SBATCH -J TensoreQTL

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 24:40:00

# mälu
#SBATCH --mem=1G

#siia alla on vaja kirjutada protsessid


module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4




NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/CEDAR/CEDAR_naive_platelet --vcf /gpfs/space/home/a82371/CEDAR/CEDAR.filtered --covariates /gpfs/space/home/a82371/CEDAR/CEDAR_microarray_platelet.covariates.txt --expression_file /gpfs/space/home/a82371/CEDAR/CEDAR_platelet_gene_expression.tsv --study CEDAR_naive_platelet --sample_genotype_ids /gpfs/space/home/a82371/CEDAR/CEDAR_naive_platelet.txt --only_autosomal_chr false -profile tartu_hpc -resume
