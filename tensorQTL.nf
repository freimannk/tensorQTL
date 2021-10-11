#!/usr/bin/env nextflow


params.covariates = ''
params.expression_file = ''  
params.outputpath = ''
params.vcf='' 
params.study=''
params.pvalue=1
params.sample_genotype_ids=''
params.only_autosomal_chr=true

// for hg38
covariates_ch = Channel.fromPath(params.covariates).collect()
expression_ch = Channel.fromPath(params.expression_file).collect()
sample_genotype_ids_ch = Channel.fromPath(params.sample_genotype_ids).collect()
sample_genotype_ids_for_filtering_ch = Channel.fromPath(params.sample_genotype_ids).collect()

variant_path = ''

if( (params.only_autosomal_chr == true) ){
  Channel
       variant_path = "$baseDir/data/variants_regions.tsv"
} else{
       variant_path = "$baseDir/data/variants_regions_with_chrX.tsv"
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
   

        """
            python3 $baseDir/scripts/generate_ge.py -g ${ge_file} -s ${ids_file} -t $baseDir/data/Homo_sapiens_GRCh38_96_genes_TSS.tsv -n ${params.study}


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

process MakeBFiles {
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'


    input:
    file(vcf) from variant_vcf_ch.flatten()


    output:
    tuple file("${vcf.simpleName}.bed"), file("${vcf.simpleName}.bim"), file("${vcf.simpleName}.fam") into bFiles_ch


    script:
        """
         plink2 --make-bed --output-chr chrM --vcf ${vcf} --out ${vcf.simpleName} --const-fid


        """

}


process TensorQTL {
    publishDir "${params.outputpath}/${params.study}", mode: 'copy', overwrite: true
    container= 'quay.io/eqtlcatalogue/tensorqtl:v21.10.1'
    input:
    set file(bed), file(bim), file(fam), file(expressionFile) from bFiles_ch.combine(tabixed_formated_ge_bed_ch.flatten())
    file covariatesFile from covariates_ch


    output:
    file '*.parquet' 


    script:
  
        """
        python3 -m tensorqtl ${bed.baseName} ${expressionFile} ${params.study}_${bed.simpleName} --covariates ${covariatesFile} --mode trans --pval_threshold ${params.pvalue} --batch_size 10000 

        """
}

   

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}


