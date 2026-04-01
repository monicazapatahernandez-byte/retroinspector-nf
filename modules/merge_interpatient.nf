process MERGE_INTERPATIENT {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcfs
    path csis

    output:
    tuple path("${params.all_prefix}.merged.vcf.gz"), path("${params.all_prefix}.merged.vcf.gz.csi")

    script:
    """
    python3 ${projectDir}/bin/merge.py \
        -samples ${vcfs.collect { it.simpleName }.join(' ')} \
        -vcf ${vcfs} \
        -o ${params.all_prefix}.merged.vcf \
        -d ${params.dist_inter}

    bcftools sort ${params.all_prefix}.merged.vcf \
        -O z -o ${params.all_prefix}.merged.vcf.gz

    bcftools index ${params.all_prefix}.merged.vcf.gz
    """
}
