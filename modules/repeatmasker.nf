process REPEATMASKER {
    conda params.genome == "t2t" \
        ? "${projectDir}/env_repeatmasker_t2t.yaml" \
        : "${projectDir}/env.yaml"

    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    tuple path(vcf), path(csi)
    path rm_lib

    output:
    path "${params.all_prefix}.merged.rm.bed"

    script:
    def lib_arg = rm_lib.name != "NO_RM_LIB" ? "${rm_lib}" : ""
    """
    mkdir -p repeatmasker

    python3 ${projectDir}/bin/checkInsertionsMultiSample.py \
        ${vcf} \
        repeatmasker \
        ${params.species} \
        ${task.cpus} \
        ${lib_arg}

    # Buscar elementos SVA F1
    python3 ${projectDir}/bin/svaf.py \
        repeatmasker/${params.all_prefix}.merged.rm.unproc.bed \
        ${projectDir}/data/sva_types.fa \
        ${params.all_prefix}.merged.rm.bed
    """
}
