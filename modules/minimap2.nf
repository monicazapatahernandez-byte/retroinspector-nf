process MINIMAP2_ALIGN {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/alns", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    minimap2 \
        -x map-ont \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
        --MD -a \
        -t ${task.cpus} \
        ${reference} ${fastq} | \
    samtools sort -@ 16 -O BAM -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}
