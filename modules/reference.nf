process GET_REFERENCE_REPEATS {
    conda "${projectDir}/env.yaml"
    executor 'local'   
    storeDir "${params.outdir}/data"

    output:
    tuple path("repeatsReferenceTE.bed.gz"), path("repeatsReferenceTE.bed.gz.csi")

    script:
    """
    wget -O hg38.fa.out.gz \
        http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz

    zcat hg38.fa.out.gz | tail -n +4 | \
        grep -E "(S|L)INE|Retroposon|LTR|DNA" | \
        awk -v'OFS=\t' '{print \$5,\$6,\$7,\$9,\$10,\$11}' > repeatsReferenceTE.bed

    bgzip -k repeatsReferenceTE.bed
    tabix -p bed --csi repeatsReferenceTE.bed.gz
    """
}

process GET_REFERENCE_GENOME {
    conda "${projectDir}/env.yaml"
    storeDir "${params.outdir}/data"

    output:
    path "GRCh38.fa"

    script:
    """
    wget -O ref.fna.gz \
        ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    gunzip -c ref.fna.gz > GRCh38.fa
    """
}

process GET_REFERENCE_T2T {
    conda "${projectDir}/env.yaml"
    storeDir "${params.outdir}/data"

    output:
    path "chm13v2.0.fa"

    script:
    """
    wget -O chm13v2.0.fa.gz \
        https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

    gunzip -c chm13v2.0.fa.gz > chm13v2.0.fa
    """
}

process GET_REFERENCE_REPEATS_T2T {
    conda "${projectDir}/env.yaml"
    storeDir "${params.outdir}/data"

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

    RM_BIN=\$(dirname \$(which RepeatMasker))
    RM_SHARE=\${RM_BIN}/../share/RepeatMasker

    cp -r \${RM_SHARE}/Libraries/ Libraries/

    \${RM_BIN}/python3 \${RM_SHARE}/famdb.py \
        -i Libraries/famdb \
        append humanAutoXYape.embl \
        --name 'T2T_CHM13_repeats'
    """
}




process PREPARE_TE_NOOVERLAP {
    conda "${projectDir}/env.yaml"
    storeDir "${params.outdir}/data"

    input:
    tuple path(repeats_bed), path(repeats_csi)

    output:
    tuple path("repeatsReferenceTENoOverlap.bed.gz"), path("repeatsReferenceTENoOverlap.bed.gz.csi")

    script:
    """
    bedtools merge -c 4,5,6 -o collapse,collapse,collapse \
        -i ${repeats_bed} > repeatsReferenceTENoOverlap.bed

    bgzip -k repeatsReferenceTENoOverlap.bed
    tabix -p bed --csi repeatsReferenceTENoOverlap.bed.gz
    """
}

process PREPARE_TE_NOOVERLAP_T2T {
    conda "${projectDir}/env.yaml"
    storeDir "${params.outdir}/data"

    input:
    tuple path(repeats_bed), path(repeats_csi)

    output:
    tuple path("repeatsReferenceTENoOverlap_t2t.bed.gz"), path("repeatsReferenceTENoOverlap_t2t.bed.gz.csi")

    script:
    """
    bedtools merge -c 4,5,6 -o collapse,collapse,collapse \
        -i ${repeats_bed} > repeatsReferenceTENoOverlap_t2t.bed

    bgzip -k repeatsReferenceTENoOverlap_t2t.bed
    tabix -p bed --csi repeatsReferenceTENoOverlap_t2t.bed.gz
    """
}


process GET_DFAM_PARTITION7 {
    conda "${projectDir}/env_repeatmasker_t2t.yaml"
    executor 'local'
    storeDir "${params.outdir}/dfam"

    output:
    path "dfam39_full*.h5"

    script:
    // Dfam 3.9 stores families partitioned by taxonomy. Partition 0 (root)
    // holds the taxonomy backbone and is mandatory for any famdb.py lookup.
    // Partition 7 contains Mammalia, which is where Homo sapiens resides.
    // Both are required for RepeatMasker to annotate human retrotransposons.
    // Format is FamDB 2.x — compatible with the famdb.py bundled in
    // RepeatMasker 4.1.9. Dfam 4.0 / FamDB 3.0 would require a separate
    // FamDB standalone package and a rebuilt image.
    """
    wget https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.0.h5.gz
    wget https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.7.h5.gz

    gzip -t dfam39_full.0.h5.gz
    gzip -t dfam39_full.7.h5.gz

    gunzip -f dfam39_full.0.h5.gz
    gunzip -f dfam39_full.7.h5.gz
    """
}
