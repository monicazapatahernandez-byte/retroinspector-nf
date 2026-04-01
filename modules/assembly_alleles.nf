process ASSEMBLY_ALLELES {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants/polished", mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf), path(csi), path(bam), path(bai)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.polished.vcf.gz"), path("${sample_id}.${caller}.polished.vcf.gz.csi")

    script:
    """
    python3 ${projectDir}/bin/getGoodAlts.py \
        ${vcf} \
        ${bam} \
        ${sample_id}.${caller}.polished.vcf \
        ${task.cpus}

    bcftools sort ${sample_id}.${caller}.polished.vcf \
        -O z -o ${sample_id}.${caller}.polished.vcf.gz

    bcftools index ${sample_id}.${caller}.polished.vcf.gz
    """
}
