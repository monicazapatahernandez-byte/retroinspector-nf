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

process GET_REFERENCE_REPEATS_T2T {
    conda "${projectDir}/env.yaml"
    publishDir "${params.outdir}/data", mode: 'copy'
    output:
    tuple path("repeatsReferenceTE_t2t.bed.gz"), path("repeatsReferenceTE_t2t.bed.gz.csi")
    script:
    """
    wget -O chm13_rm_raw.bed \
        https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
    sort -k1,1 -k2,2n chm13_rm_raw.bed | bgzip > repeatsReferenceTE_t2t.bed.gz
    tabix --csi -p bed repeatsReferenceTE_t2t.bed.gz
    """
}

process GET_ANNOTATION_T2T {
    conda "${projectDir}/env.yaml"
    publishDir "${params.outdir}/data", mode: 'copy'
    output:
    path "hs1.ncbiRefSeq.gtf"
    script:
    """
    wget -O hs1.ncbiRefSeq.gtf.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz
    gunzip hs1.ncbiRefSeq.gtf.gz
    """
}

process GET_REPEATMASKER_LIB_T2T {
    conda "${projectDir}/env_repeatmasker_t2t.yaml"
    storeDir "${params.outdir}/data/rm_lib_t2t"
    output:
    path "Libraries"
    script:
    """
    wget -O humanAutoXYape.embl \
        https://raw.githubusercontent.com/jessicaStorer88/RepeatMasker_library_CHM13/main/humanAutoXYape.embl
    RM_SHARE=\$(dirname \$(which RepeatMasker))/../share/RepeatMasker
    cp -r \${RM_SHARE}/Libraries/ Libraries/
    python3 \${RM_SHARE}/famdb.py \
        -i Libraries/famdb \
        append humanAutoXYape.embl \
        --name 'T2T_CHM13_repeats'
    """
}
