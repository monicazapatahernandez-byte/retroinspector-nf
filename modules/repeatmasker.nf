process REPEATMASKER {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    tuple path(vcf), path(csi)

    output:
    path "${params.all_prefix}.merged.rm.bed"

    script:
    """
    mkdir -p repeatmasker

    python3 ${projectDir}/bin/checkInsertionsMultiSample.py \
        ${vcf} \
        repeatmasker \
        ${params.species} \
        ${task.cpus}

    # Buscar elementos SVA F1
    python3 ${projectDir}/bin/svaf.py \
        repeatmasker/${params.all_prefix}.merged.rm.unproc.bed \
        ${projectDir}/data/sva_types.fa \
        ${params.all_prefix}.merged.rm.bed
    """
}
