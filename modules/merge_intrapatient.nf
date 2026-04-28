process MERGE_INTRAPATIENT {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    input:
    tuple val(sample_id), path(vcfs), path(csis)

    output:
    tuple val(sample_id), path("${sample_id}.merged.both.ungenotyped.vcf")

    script:
    """
    python3 ${projectDir}/bin/merge.py \
        -samples ${params.callers.tokenize(',').join(' ')}  \
        -vcf ${vcfs} \
        -o ${sample_id}.merged.both.ungenotyped.vcf \
        -d ${params.dist_intra}
    """
}
