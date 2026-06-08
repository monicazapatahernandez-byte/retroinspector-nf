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


    def dfam_list = dfam_files instanceof List
        ? dfam_files
        : (dfam_files.name != "NO_DFAM7" ? [dfam_files] : [])

    def libdir_setup = dfam_list ? """mkdir -p rm_libdir/famdb
   for entry in ${lib_arg}/*; do
       name=\$(basename "\$entry")
       [ "\$name" = "famdb" ] && continue
       ln -sf \$(realpath "\$entry") rm_libdir/
   done
   for f in ${dfam_list.join(' ')}; do
       ln -sf \$(realpath \$f) rm_libdir/famdb/
   done""" : ""


    def libdir_for_rm = dfam_list ? "rm_libdir" : lib_arg

    """
    export PATH="${env_dir}/bin:\${PATH}"
    ${libdir_setup}

    mkdir -p repeatmasker

    ${env_dir}/bin/python3 ${projectDir}/bin/checkInsertionsMultiSample.py \\
        ${vcf} \\
        repeatmasker \\
        ${params.species} \\
        ${task.cpus} \\
        ${libdir_for_rm}

    
    ${env_dir}/bin/python3 ${projectDir}/bin/svaf.py \\
        repeatmasker/${params.all_prefix}.merged.rm.unproc.bed \\
        ${projectDir}/data/sva_types.fa \\
        ${params.all_prefix}.merged.rm.bed
    """
}
