process ASSEMBLY_ALLELES {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants/polished", mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf), path(csi), path(bam), path(bai)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.polished.vcf.gz"), path("${sample_id}.${caller}.polished.vcf.gz.csi")

    script:
    """
    N_VARS=\$(zcat ${vcf} | grep -v "^#" | wc -l)
    echo "Variantes en ${vcf}: \$N_VARS"

    if [ "\$N_VARS" -eq 0 ]; then
        echo "AVISO: VCF vacío, generando salida vacía válida"
        zcat ${vcf} | grep "^#" > ${sample_id}.${caller}.polished.vcf
    else
        python3 ${projectDir}/bin/getGoodAlts.py \
            ${vcf} \
            ${bam} \
            ${sample_id}.${caller}.polished.vcf \
            ${task.cpus}
    fi

    bcftools sort ${sample_id}.${caller}.polished.vcf \
        -O z -o ${sample_id}.${caller}.polished.vcf.gz

    bcftools index ${sample_id}.${caller}.polished.vcf.gz
    """
}
