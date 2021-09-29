# tensorQTL
tensorQTL workflow
NXF_VER=20.10.0 ./nextflow run tensorQTL.nf --outputpath {outputpath}
                                            --vcf {vcf file} 
                                            --covariates {covariates file} 
                                            --expression_file {expression file}   
                                            --study {study name} 
                                            --sample_genotype_ids {file with sample id and matching genotype id}  
                                            --variant_ranges {file with variant position ranges to analyse in one process}
                                            --genes_tss {file to  specify the center of the cis-window (usually the TSS), with start == end-1}
                                            -resume
