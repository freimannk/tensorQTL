# tensorQTL
tensorQTL workflow
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath {outputpath} <br />
                                            --vcf {vcf file} <br />
                                            --covariates {covariates file} <br />
                                            --expression_file {expression file}  <br /> 
                                            --study {study name} <br />
                                            --sample_genotype_ids {file with sample id and matching genotype id}  <br />
                                            --variant_ranges {file with variant position ranges to analyse in one process} <br />
                                            --genes_tss {file to  specify the center of the cis-window (usually the TSS), with start == end-1} <br />
                                            -resume
