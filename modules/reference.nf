process GET_REFERENCE_REPEATS {
    conda "${projectDir}/env.yaml"
    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    tuple path("repeatsReferenceTE.bed.gz"), path("repeatsReferenceTE.bed.gz.csi")

    script:
    """
    cp ${projectDir}/data/repeatsReferenceTE.bed.gz repeatsReferenceTE.bed.gz
    cp ${projectDir}/data/repeatsReferenceTE.bed.gz.csi repeatsReferenceTE.bed.gz.csi
    """
}

process GET_REFERENCE_GENOME {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    path "GRCh38.fa"

    script:
    """
    wget -O ref.fna.gz \
        ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    gunzip ref.fna.gz > GRCh38.fa
    """
}

process GET_REFERENCE_T2T {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    path "chm13v2.0.fa"

    script:
    """
    wget -O chm13v2.0.fa.gz \
        https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
    gunzip chm13v2.0.fa.gz
    """
}
