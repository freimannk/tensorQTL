#!/usr/bin/env nextflow

nextflow.enable.dsl = 1

params.covariates = ''
params.expression_file = ''  
params.outputpath = ''
params.vcf_prefix_path='' 
params.dataset=''  
params.pvalue=1
params.sample_genotype_ids=''
params.only_autosomal_chr=true
params.median_tpm_filtration_file=''
params.vcf_genotype_field='DS'
params.maf_filter=0.05
params.hwe=1e-5
params.testing = false


 

covariates_ch = Channel.fromPath(params.covariates).collect()
expression_ch = Channel.fromPath(params.expression_file).collect()
sample_genotype_ids_ch = Channel.fromPath(params.sample_genotype_ids).collect()
sample_genotype_ids_for_filtering_ch = Channel.fromPath(params.sample_genotype_ids).collect()
sample_genotype_ids_for_cov_filtering_ch = Channel.fromPath(params.sample_genotype_ids).collect()


variant_path = ''
if(params.testing == true){
    variant_path= "$baseDir/data/testing_data/variants_regions_testing.tsv" 
} 
else if(params.only_autosomal_chr == true){
    variant_path="$baseDir/data/variants_regions.tsv" 
} else{
    variant_path = "$baseDir/data/variants_regions_with_chrX.tsv"
}

Channel
       .fromPath(variant_path)
       .splitText().map{it -> it.trim()}
       .set{variant_ranges_ch} 


Channel
    .from(params.vcf_prefix_path)
    .map { dataset -> [file("${dataset}.vcf.gz"), file("${dataset}.vcf.gz.csi")]}
    .set { vcf_ch }


process FilterSamplesFromVCF {
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    set file(huge_vcf), file(indexed) from vcf_ch
    file ids_file from sample_genotype_ids_for_filtering_ch


    output:
    tuple file("${huge_vcf.simpleName}_tagged_filtered.vcf.gz"), file("${huge_vcf.simpleName}_tagged_filtered.vcf.gz.csi") into filtered_vcf_ch     


    script:
        """
        awk '(NR>1)''{print \$2}' ${ids_file} > ${huge_vcf.simpleName}.txt
        bcftools view -S ${huge_vcf.simpleName}.txt ${huge_vcf} -Oz -o ${huge_vcf.simpleName}___samples_filtered.vcf.gz
        bcftools +fill-tags ${huge_vcf.simpleName}___samples_filtered.vcf.gz -Oz -o ${huge_vcf.simpleName}_tagged.vcf.gz
        bcftools filter -i 'INFO/HWE > ${params.hwe} & MAF[0] > ${params.maf_filter}' ${huge_vcf.simpleName}_tagged.vcf.gz -Oz -o ${huge_vcf.simpleName}_tagged_filtered.vcf.gz
        bcftools index ${huge_vcf.simpleName}_tagged_filtered.vcf.gz
        """

}


process PrepateGeneExpressionFile {
    container = 'quay.io/fhcrc-microbiome/python-pandas:4a6179f'

    input:
    file ge_file from expression_ch
    file ids_file from sample_genotype_ids_ch



    output:
    file '*_tensorQTL_wf_ge.bed' into formated_ge_bed_ch 



    script:
    median_tpm_filtration_file = params.median_tpm_filtration_file ? "-f ${params.median_tpm_filtration_file}" : ""

   

        """
            generate_ge.py -g ${ge_file} -s ${ids_file} -t $baseDir/data/Homo_sapiens_GRCh38_96_genes_TSS.tsv -n ${params.dataset} ${median_tpm_filtration_file}


        """

}

process TabixGEBedFile {
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    file ge_file from formated_ge_bed_ch



    output:
    file '*.bed.gz' into tabixed_formated_ge_bed_ch 



    script:
   
 
        """
            bgzip ${ge_file} 
            tabix -p bed ${ge_file}.gz

        """

}


process GenerateVariantVCFFiles {
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    set file(vcf), file(indexed), val(variant_range) from filtered_vcf_ch.combine(variant_ranges_ch.flatten())

    output:
        file("chr+${chr_str}_${range_str}+.vcf.gz") into variant_vcf_ch  


    script:
    variant_range_string = variant_range.split(':')
    chr_str = variant_range_string[0]
    range_str = variant_range_string[1]
        """
         bcftools view ${vcf} -r ${variant_range} -t ${variant_range} -Oz -o chr+${chr_str}_${range_str}+.vcf.gz 
        """

}



process VcfToDosage{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    input:
    file(vcf) from variant_vcf_ch.flatten()

    output:
    file("${vcf.simpleName}.tsv.gz") into genotype_matrix_DS_ch

    script:
    if(params.vcf_genotype_field == 'DS'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nID\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz
        #Extract dosage and merge
        bcftools query -f "%CHROM\\t%POS\\t%ID[\\t%DS]\\n" ${vcf} | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.tsv.gz
        """
    } else if (params.vcf_genotype_field == 'GT'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nID\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz
        #Extract dosage and merge
        bcftools +dosage ${vcf} -- -t GT | tail -n+2 | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.tsv.gz
        """
    } 
}

    process FilterCovariates{
    container = 'quay.io/fhcrc-microbiome/python-pandas:4a6179f'

    input:
    file covariatesFile from covariates_ch
    file ids_file from sample_genotype_ids_for_cov_filtering_ch


    output:
    file("${params.dataset}_filtered_covariates.txt") into filtered_covariates_ch

    script:

        """
        filter_covariates.py -c ${covariatesFile} -n ${params.dataset} -s ${ids_file}
        """
    
}
        
    process CorrectMissingDosage{
    container = 'quay.io/fhcrc-microbiome/python-pandas:4a6179f'


    input:
    file(matrix_csv) from genotype_matrix_DS_ch.flatten()

    output:
    file("${matrix_csv.simpleName}.txt.gz") into genotype_matrix_DS_missing_removed_ch

    script:
        """
        correct_missing_DS.py -i ${matrix_csv} -n ${matrix_csv.simpleName}

        """
}



process TensorQTL {
    publishDir "${params.outputpath}/${params.dataset}", mode: 'copy', overwrite: true
    container= 'quay.io/eqtlcatalogue/tensorqtl:v21.10.1'
    input:
    set file(genotype_matrix), file(expressionFile) from genotype_matrix_DS_missing_removed_ch.combine(tabixed_formated_ge_bed_ch.flatten())
    file covariatesFile from filtered_covariates_ch


    output:
    file '*.parquet' optional true


    script:
  
        """
        tensorQTL.py -e ${expressionFile} -c ${covariatesFile} -g ${genotype_matrix} -p ${params.pvalue} -m ${params.maf_filter} -n ${params.dataset}_${genotype_matrix.simpleName}
        """
}

   

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}


