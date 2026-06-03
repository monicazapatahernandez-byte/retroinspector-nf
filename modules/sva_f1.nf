process SVA_F1 {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    path rm_bed

    output:
    path "rm_svaf1.bed"

    script:
    """
    export PATH="/opt/conda/envs/retro-base/bin:\${PATH}"

    /opt/conda/envs/retro-base/bin/python3 ${projectDir}/bin/svaf.py \
        ${rm_bed} \
        rm_svaf1.bed \
        svaf.fasta \
        svaf.blast6 \
        svaf_db \
        ${projectDir}/data/sva_types.fa
    """
}
