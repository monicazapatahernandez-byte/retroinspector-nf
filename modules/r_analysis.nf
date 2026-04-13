process R_PREPARATORY {
    conda "${projectDir}/r.yaml"

    publishDir "${params.outdir}/rds", mode: 'copy'

    input:
    path rm_bed
    tuple path(del_vcf), path(del_csi)

    output:
    path "annotation.rds"
    path "annotatedInsertionsMin3.rds"
    path "insertionsTable.rds"
    path "allIns.rds"
    path "meDeletionsMin3.rds"

    script:
    """
    Rscript ${projectDir}/bin/analysisPreparatory.R \
        ${rm_bed} \
        ${del_vcf} \
        ${params.min_read_support} \
        ${params.all_prefix}
    """
}

process R_GENOTYPING {
    conda "${projectDir}/r.yaml"

    publishDir "${params.outdir}/rds", mode: 'copy'

    input:
    path insertions_table
    path annotated_insertions

    output:
    path "genes.rds"
    path "${params.all_prefix}.me.insertions.txt"
    path "${params.all_prefix}.me.insertions.lax.txt"

    script:
    """
    Rscript ${projectDir}/bin/analysisGenotyping.R \
        ${insertions_table} \
        ${annotated_insertions} \
        ${params.all_prefix}
    """
}

process R_ENRICHMENT {
    conda "${projectDir}/r.yaml"

    publishDir "${params.outdir}/rds", mode: 'copy'

    input:
    path genes

    output:
    path "egoMF.rds"
    path "egoBP.rds"
    path "egoCC.rds"
    path "do.rds"
    path "ncg.rds"

    script:
    """
    Rscript ${projectDir}/bin/enrichment.R \
        ${genes} \
        ${params.enrichment_pval}
    """
}

process R_REPORT {
    conda "${projectDir}/r.yaml"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path insertions_table
    path all_ins
    path annotated_insertions
    path me_deletions
    path egoMF
    path egoBP
    path egoCC
    path do_rds
    path ncg
    path other_sets

    output:
    path "report.${params.all_prefix}.html"

    script:
    """
    Rscript -e "rmarkdown::render('${projectDir}/bin/report.Rmd',
        output_file='report.${params.all_prefix}.html',
        output_dir=getwd(),
        params=list(
            insertions_table='${insertions_table}',
            all_ins='${all_ins}',
            annotated_insertions='${annotated_insertions}',
            me_deletions='${me_deletions}',
            egoMF='${egoMF}',
            egoBP='${egoBP}',
            egoCC='${egoCC}',
            do_rds='${do_rds}',
            ncg='${ncg}',
            plimit=${params.enrichment_pval}
        ))"
    """

process R_COMPARE {
    conda "${projectDir}/r.yaml"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample1), val(sample2)
    path insertions_table
    path annotated_insertions

    output:
    path "${sample1}_vs_${sample2}.html"

    script:
    """
    Rscript -e "rmarkdown::render('${projectDir}/bin/comparison.Rmd',
        output_file='${sample1}_vs_${sample2}.html',
        output_dir=getwd(),
        params=list(
            insertions_table='${insertions_table}',
            annotated_insertions='${annotated_insertions}',
            sample1='${sample1}',
            sample2='${sample2}',
            samples='${params.all_prefix}'
        ))"
    """
}

process GENERATE_VCF {
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcf_body
    path vcf_body_lax
    path header

    output:
    tuple path("${params.all_prefix}.te.vcf.gz"), path("${params.all_prefix}.te.vcf.gz.csi")
    tuple path("${params.all_prefix}.te.lax.vcf.gz"), path("${params.all_prefix}.te.lax.vcf.gz.csi")

    script:
    """
    cat ${header} ${vcf_body} | \
        bcftools sort -O z -o ${params.all_prefix}.te.vcf.gz
    bcftools index ${params.all_prefix}.te.vcf.gz

    cat ${header} ${vcf_body_lax} | \
        bcftools sort -O z -o ${params.all_prefix}.te.lax.vcf.gz
    bcftools index ${params.all_prefix}.te.lax.vcf.gz
    """
}
