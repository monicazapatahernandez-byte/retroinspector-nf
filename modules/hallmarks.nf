process HALLMARKS {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/rds", mode: 'copy'

    input:
    path hallmarks_input
    path reference

    output:
    path "hallmarks.tsv"

    script:
    """
    export PATH="/opt/conda/envs/retro-base/bin:\${PATH}"

    /opt/conda/envs/retro-base/bin/samtools faidx ${reference}

    /opt/conda/envs/retro-base/bin/python3 ${projectDir}/bin/hallmarks.py \
        --input     ${hallmarks_input} \
        --reference ${reference} \
        --output    hallmarks.tsv
    """
}
