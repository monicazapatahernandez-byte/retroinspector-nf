process MERGE_INTRAPATIENT {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    input:
    tuple val(sample_id), path(vcfs), path(csis)

    output:
    tuple val(sample_id), path("${sample_id}.merged.both.ungenotyped.vcf")

    script:
    """
    export PATH="/opt/conda/envs/retro-base/bin:\${PATH}"

    CUTESV_VCF=\$(ls ${sample_id}.cutesv.polished.vcf.gz)
    SNIFFLES2_VCF=\$(ls ${sample_id}.sniffles2.polished.vcf.gz)

    /opt/conda/envs/retro-base/bin/python3 ${projectDir}/bin/merge.py \
        -samples cuteSV sniffles2 \
        -vcf \${CUTESV_VCF} \${SNIFFLES2_VCF} \
        -o ${sample_id}.merged.both.ungenotyped.vcf \
        -d ${params.dist_intra}
    """
}
