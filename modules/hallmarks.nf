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
    samtools faidx ${reference}

    \${CONDA_PREFIX}/bin/python3 ${projectDir}/bin/hallmarks.py \
        --input     ${hallmarks_input} \
        --reference ${reference} \
        --output    hallmarks.tsv
    """
}
