#!/bin/bash

#Töö nimi
#SBATCH -J TEST

#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#Töö kestus
#SBATCH -t 03:40:00

# mälu
#SBATCH --mem=5G

module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load nextflow
module load squashfs/4.4


NXF_VER=20.10.0 nextflow run tensorQTL.nf -profile tartu_hpc,testing_tartu_hpc -resume
