process SURVIVOR_INTRAPATIENT {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    input:
    tuple val(sample_id), path(vcfs), path(csis)

    output:
    tuple val(sample_id), path("${sample_id}.merged.survivor.ungenotyped.vcf")

    script:
    """
    surpyvor merge \
        --variants ${vcfs} \
        -d ${params.survivor_intra} \
        -c 1 -l 1 \
        -o ${sample_id}.merged.survivor.ungenotyped.vcf
    """
}

process MOSDEPTH_DEL {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    input:
    tuple val(sample_id), path(vcf), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.del.regions.bed.gz"), path(vcf)

    script:
    """
    python3 ${projectDir}/bin/vcfToBedForMosdepth.py \
        --svtype DEL \
        ${vcf} ${sample_id}.del.bed

    mosdepth \
        -t ${task.cpus} \
        -Q 20 -n \
        -b ${sample_id}.del.bed \
        ${sample_id}.del \
        ${bam}
    """
}

process GENOTYPE_DEL {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants/survivor", mode: 'copy'

    input:
    tuple val(sample_id), path(coverage), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.merged.survivor.vcf.gz"), path("${sample_id}.merged.survivor.vcf.gz.csi")

    script:
    """
    python3 ${projectDir}/bin/genotype.py \
        --svtype DEL \
        ${coverage} \
        ${vcf} \
        ${sample_id}.merged.survivor.vcf

    bcftools sort -O z \
        -o ${sample_id}.merged.survivor.vcf.gz \
        ${sample_id}.merged.survivor.vcf

    bcftools index ${sample_id}.merged.survivor.vcf.gz
    """
}

process SURVIVOR_INTERPATIENT {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants/survivor", mode: 'copy'

    input:
    path vcfs
    path csis

    output:
    tuple path("${params.all_prefix}.merged.survivor.vcf.gz"), path("${params.all_prefix}.merged.survivor.vcf.gz.csi")

    script:
    """
    for vcf in ${vcfs}; do
        bcftools view -O v -o \${vcf%.gz} \$vcf
        echo \${vcf%.gz} >> survivor_list.txt
    done

    SURVIVOR merge survivor_list.txt \
        ${params.survivor_inter} 1 1 -1 -1 1 \
        ${params.all_prefix}.merged.survivor.vcf

    cat ${params.all_prefix}.merged.survivor.vcf | \
        python3 ${projectDir}/bin/fixVCF.py | \
        bcftools sort -O z \
        -o ${params.all_prefix}.merged.survivor.vcf.gz

    bcftools index ${params.all_prefix}.merged.survivor.vcf.gz
    """
}

process GET_DELETIONS {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple path(vcf), path(csi)
    tuple path(bed), path(bed_csi)

    output:
    tuple path("${params.all_prefix}.me.deletions.vcf.gz"), path("${params.all_prefix}.me.deletions.vcf.gz.csi")

    script:
    """
    python3 ${projectDir}/bin/getMeDeletions.py \
        ${vcf} \
        ${bed} \
        ${params.all_prefix}.me.deletions.temp.vcf.gz

    bcftools sort \
        ${params.all_prefix}.me.deletions.temp.vcf.gz \
        -O z -o ${params.all_prefix}.me.deletions.vcf.gz

    bcftools index ${params.all_prefix}.me.deletions.vcf.gz
    """
}
