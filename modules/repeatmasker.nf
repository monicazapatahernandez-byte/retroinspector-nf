process REPEATMASKER {
    conda params.genome == "t2t" \
        ? "${projectDir}/env_repeatmasker_t2t.yaml" \
        : "${projectDir}/env.yaml"

    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    tuple path(vcf), path(csi)
    path rm_lib
    path dfam7

    output:
    path "${params.all_prefix}.merged.rm.bed"

    script:
    def lib_arg   = rm_lib.name  != "NO_RM_LIB"  ? "${rm_lib}"  : ""
    def famdb_dir = "\${CONDA_PREFIX}/share/RepeatMasker/Libraries/famdb"
    def dfam_cp = dfam7.name != "NO_DFAM7"
    ? """[ -f ${famdb_dir}/${dfam7.name} ] || ln -sf \$(realpath ${dfam7}) ${famdb_dir}/
       ln -sf \$(realpath ${dfam7}) ${lib_arg}/famdb/"""
    : ""
    """
    export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    ${dfam_cp}

    mkdir -p repeatmasker

    \${CONDA_PREFIX}/bin/python3 ${projectDir}/bin/checkInsertionsMultiSample.py \
        ${vcf} \
        repeatmasker \
        ${params.species} \
        ${task.cpus} \
        ${lib_arg}

    # Buscar elementos SVA F1
    \${CONDA_PREFIX}/bin/python3 ${projectDir}/bin/svaf.py \
        repeatmasker/${params.all_prefix}.merged.rm.unproc.bed \
        ${projectDir}/data/sva_types.fa \
        ${params.all_prefix}.merged.rm.bed
    """
}
