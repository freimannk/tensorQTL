# tensorQTL workflow
 <br />
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf <br /> --outputpath {outputpath} <br />
                                            --vcf {prefix path to vcf file and its indexed file} <br />
                                            --covariates {covariates file} <br />
                                            --expression_file {expression file}  <br /> 
                                            --study {study name} <br />
                                            --sample_genotype_ids {file with sample id and matching genotype id, with header (/example_files)}  <br />
                                            --variant_ranges {file with variant position ranges to analyse in one process} <br />
                                            --genes_tss {file to  specify the center of the cis-window (usually the TSS), with start == end-1} <br />
                                            --pvalue {to filter results by p-value, default: 1}  <br />
                                            -resume  <br />

#### genes_tss file and variant_ranges file are in the data/ directory.

 <br />
 
#### TODOs: 
<br/>

                - francois4/tensorqtl:latest docker img too big (fails sometimes)? <br/>
                - biocontainers/bcftools:v1.9-1-deb_cv1 img do not contain ps which is needed for NF stats collection for trave file
                

