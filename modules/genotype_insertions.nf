process MOSDEPTH_INS {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    input:
    tuple val(sample_id), path(vcf), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.ins.regions.bed.gz"), path(vcf)

    script:
    """
    python3 ${projectDir}/bin/vcfToBedForMosdepth.py \
        ${vcf} ${sample_id}.ins.bed

    mosdepth \
        -t ${task.cpus} \
        -Q 20 -n \
        -b ${sample_id}.ins.bed \
        ${sample_id}.ins \
        ${bam}
    """
}

process GENOTYPE_INS {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(coverage), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.merged.both.vcf.gz"), path("${sample_id}.merged.both.vcf.gz.csi")

    script:
    """
    python3 ${projectDir}/bin/genotype.py \
        ${coverage} \
        ${vcf} \
        ${sample_id}.merged.both.vcf

    bcftools sort -O z -o ${sample_id}.merged.both.vcf.gz \
        ${sample_id}.merged.both.vcf

    bcftools index ${sample_id}.merged.both.vcf.gz
    """
}
