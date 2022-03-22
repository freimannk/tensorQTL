# tensorQTL workflow
 <br />
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf <br /> --outputpath {outputpath} <br />
                                             --vcf {prefix path to vcf file and its indexed file} <br />
                                             --covariates {covariates file} <br />
                                             --expression_file {expression file}  <br /> 
                                             --study {study name} <br />
                                             --sample_genotype_ids {file with sample id and matching genotype id, with header (example is in the /example_files)}  <br />
                                             --pvalue {to filter results by p-value, default: 1}  <br />
                                             --only_autosomal_chr {boolean, if false includes regions from X chr, default: true}  <br />
                                             --median_tpm_filtration_file {median tmp file to filter genes, needed to speed up analysis, otherwise includes all genes}  <br />
                                             -profile tartu_hpc  <br />
                                             -resume  <br /> 
 

