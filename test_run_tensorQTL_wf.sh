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




NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath /gpfs/space/home/a82371/TensorQTL/blueprint_naive_Tcells/tensorQTLResults_docker --vcf /gpfs/space/home/a82371/eQTL/genotypes/merged_indexed_gen_files/all_blueprint_schmiedel_bossini-castillo_gencord_rm_missing_gt_corrected_ids_exclude_multiallelic --covariates /gpfs/space/home/a82371/tensorqtl/data/covariates/BLUEPRINT_PE.T-cell_ge.covariates.txt --expression_file /gpfs/space/home/a82371/TensorQTL/rna_seq_ge_qcnorm/naive_T_BLUE_SCHM_BOS_GENC.naive_T_cell.tsv --study blueprint_naive_Tcells --sample_genotype_ids /gpfs/space/home/a82371/TensorQTL/blueprint_naive_Tcells/blueprint_naive_Tcells_ids.txt --variant_ranges /gpfs/space/home/a82371/TensorQTL/schmiedel_naive_Tcells/testing_variants_gen_position_range.tsv --genes_tss /gpfs/space/home/a82371/tensorQTL_workflow/tensorQTL/data/Homo_sapiens_GRCh38_96_genes_TSS.tsv -profile tartu_hpc -resume
