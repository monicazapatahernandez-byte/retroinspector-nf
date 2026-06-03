process REPEATMASKER {
    conda params.genome == "t2t" \
        ? "${projectDir}/env_repeatmasker_t2t.yaml" \
        : "${projectDir}/env.yaml"

    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    tuple path(vcf), path(csi)
    path rm_lib
    path dfam_files

    output:
    path "${params.all_prefix}.merged.rm.bed"

    script:
    def lib_arg = rm_lib.name != "NO_RM_LIB" ? "${rm_lib}" : ""
    def env_dir = params.genome == "t2t" \
        ? "/opt/conda/envs/retro-rm-t2t" \
        : "/opt/conda/envs/retro-base"

    def famdb_dir = "${env_dir}/share/RepeatMasker/Libraries/famdb"

    def dfam_cp = dfam_files instanceof List
        ? """for f in ${dfam_files.join(' ')}; do
        [ -f ${famdb_dir}/\$(basename \$f) ] || ln -sf \$(realpath \$f) ${famdb_dir}/
        ln -sf \$(realpath \$f) ${lib_arg}/famdb/
   done"""
        : dfam_files.name != "NO_DFAM7"
        ? """[ -f ${famdb_dir}/${dfam_files.name} ] || ln -sf \$(realpath ${dfam_files}) ${famdb_dir}/
   ln -sf \$(realpath ${dfam_files}) ${lib_arg}/famdb/"""
        : ""

    """
    export PATH="${env_dir}/bin:\${PATH}"

    ${dfam_cp}

    mkdir -p repeatmasker

    ${env_dir}/bin/python3 ${projectDir}/bin/checkInsertionsMultiSample.py \
        ${vcf} \
        repeatmasker \
        ${params.species} \
        ${task.cpus} \
        ${lib_arg}

    # Buscar elementos SVA F1
    ${env_dir}/bin/python3 ${projectDir}/bin/svaf.py \
        repeatmasker/${params.all_prefix}.merged.rm.unproc.bed \
        ${projectDir}/data/sva_types.fa \
        ${params.all_prefix}.merged.rm.bed
    """
}
