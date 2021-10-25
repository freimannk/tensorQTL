#!/usr/bin/env nextflow


params.covariates = ''
params.expression_file = ''  
params.outputpath = ''
params.vcf='' 
params.study=''
params.pvalue=1
params.sample_genotype_ids=''
params.only_autosomal_chr=true
params.median_tpm_filtration_file=''
params.vcf_genotype_field='DS'
params.filter_cis_window=''


// for hg38
covariates_ch = Channel.fromPath(params.covariates).collect()
expression_ch = Channel.fromPath(params.expression_file).collect()
sample_genotype_ids_ch = Channel.fromPath(params.sample_genotype_ids).collect()
sample_genotype_ids_for_filtering_ch = Channel.fromPath(params.sample_genotype_ids).collect()

variant_path = ''

if( (params.only_autosomal_chr == true) ){
       variant_path = "$baseDir/data/variant_regions_containing_500_variants.tsv"
} else{
       variant_path = "$baseDir/data/variants_regions_w_chrX.tsv"
}

Channel
       .fromPath(variant_path)
       .splitText().map{it -> it.trim()}
       .set{variant_ranges_ch} 


Channel
    .from(params.vcf)
    .map { study -> [file("${study}.vcf.gz"), file("${study}.vcf.gz.csi")]}
    .set { vcf_ch }


process FilterSamplesFromVCF {
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    set file(huge_vcf), file(indexed) from vcf_ch
    file ids_file from sample_genotype_ids_for_filtering_ch


    output:
    tuple file("${huge_vcf.simpleName}___filtered.vcf.gz"), file("${huge_vcf.simpleName}___filtered.vcf.gz.csi") into filtered_vcf_ch

    script:
        """
        awk '(NR>1)''{print \$2}' ${ids_file} > ${huge_vcf.simpleName}.txt
        bcftools view -S ${huge_vcf.simpleName}.txt ${huge_vcf} -Oz -o ${huge_vcf.simpleName}___filtered.vcf.gz
        bcftools index ${huge_vcf.simpleName}___filtered.vcf.gz
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
            python3 $baseDir/scripts/generate_ge.py -g ${ge_file} -s ${ids_file} -t $baseDir/data/Homo_sapiens_GRCh38_96_genes_TSS.tsv -n ${params.study} ${median_tpm_filtration_file}


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
    file("${vcf.simpleName}.tsv.gz") into genotype_matrix_ch

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



process TensorQTL {
    publishDir "${params.outputpath}/${params.study}", mode: 'copy', overwrite: true
    container= 'quay.io/eqtlcatalogue/tensorqtl:v21.10.1'
    input:
    set file(genotype_matrix), file(expressionFile) from genotype_matrix_ch.combine(tabixed_formated_ge_bed_ch.flatten())
    file covariatesFile from covariates_ch


    output:
    file '*.parquet' optional true


    script:
        filter_cis_window = params.filter_cis_window ? "-f ${params.filter_cis_window}" : ""

  
        """
        python3 $baseDir/scripts/tensorQTL.py -e ${expressionFile} -c ${covariatesFile} -g ${genotype_matrix} -p ${params.pvalue} -n ${params.study}_${genotype_matrix.simpleName} ${filter_cis_window}
        """
}

   

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}


