# tensorQTL workflow
**Portable QTL mapper**

tensorQTL workflow uses GPU-enabled tensorqtl QTL mapper: [Taylor-Weiner, Aguet, et al., *Genome Biol.* 20:228, 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1836-7) 


## Usage examples

###  For testing

```
nextflow run tensorQTL.nf -profile tartu_hpc,testing_tartu_hpc -resume
```

Pipeline should create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
testing_results/GEUVADIS_test_samples         # Finished testing results (containing one 204M sized .parquet file)
```

###  Running the analysis

```
nextflow run tensorQTL.nf -profile tartu_hpc -resume\
 --dataset GEUVADIS\
 --vcf_prefix_path ~/tensorQTL_workflow/tensorQTL/data/testing_data/GEUVADIS_chr1_212626194_212851227_testing_samples\
 --covariates ~/tensorQTL_workflow/tensorQTL/data/testing_data/GEUVADIS_LCL_cells.covariates_added_sex_test_samples.txt\
 --expression_file ~/tensorQTL_workflow/tensorQTL/data/testing_data/GEUVADIS.LCL_naive_test_samples.tsv\
 --sample_genotype_ids ~/tensorQTL_workflow/tensorQTL/data/testing_data/GEUVADIS_LCL_study_ids_test_samples.txt\
 --only_autosomal_chr false\
 --outputpath ~/results\
 --median_tpm_filtration_file ~/tensorQTL_workflow/tensorQTL/data/LCL_naive_median_tpm_filtered_genes.tsv.gz
```

Optional arguments/info

* --vcf_prefix_path : path to vcf file and its indexed file, example: `~/BLUEPRINT/BLUEPRINT.MAF001`, /BLUEPRINT folder has to contain files: `BLUEPRINT.MAF001.vcf.gz` and `BLUEPRINT.MAF001.vcf.gz.csi`
* --covariates : covariates file, from qtlmap workflow
* --vcf_genotype_field : `DS` or `GT`, default= `DS`
* --maf_filter : MAF filter, default = 0.05
* --sample_genotype_ids : file with sample id and matching genotype id, with header
* --pvalue : to filter out results by p-value, default: 1
* --only_autosomal_chr : boolean, if false includes regions from X chr, default: true
* --median_tpm_filtration_file : median tpm file to filter genes, needed to speed up analysis, otherwise includes all genes. Examples are: data/LCL_naive_median_tpm_filtered_genes.tsv.gz and data/T-cell_median_tpm_filtered_genes.tsv.gz, from qcnorm workflow
