process CUTESV {
    tag "$sample_id"
    conda "${projectDir}/env.yaml"

    publishDir "${params.outdir}/variants/cutesv", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.cutesv.vcf.gz"), path("${sample_id}.cutesv.vcf.gz.csi")

        script:
    def fixScript = "${projectDir}/bin/fixVCF.py"
    """
    BASE_ENV="/opt/conda/envs/retro-base"
    export PATH="\${BASE_ENV}/bin:\${PATH}"

       mkdir -p tmp_cutesv

    \${BASE_ENV}/bin/cuteSV \
        -t ${task.cpus} \
        --max_cluster_bias_INS 100 \
        -s 2 -L -1 \
        --report_readid \
        --sample ${sample_id} \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        ${bam} ${reference} \
        ${sample_id}.cutesv.withseq.vcf tmp_cutesv

    # Filter INS y BND, then read support
    cat ${sample_id}.cutesv.withseq.vcf | \
    gawk -v 'OFS=\t' \
    '{if (substr(\$0,1,1)=="#") {print} \
    else { match(\$8,/SVTYPE=([^;]+)/,a); \
      if (a[1]=="INS"||a[1]=="BND") {print} \
      else {\$5="<"a[1]">";\$4="N";print} }}' | \
    gawk -v 'OFS=\t' \
    '{if (substr(\$0,1,1)=="#") {print} \
    else { match(\$8,/RE=([0-9]+)/,a); \
      if (a[1]+0 >= ${params.min_read_support}) {print} }}' | \
    \${BASE_ENV}/bin/python3 ${fixScript} | \
    bcftools sort -O z -o ${sample_id}.cutesv.vcf.gz

    bcftools index ${sample_id}.cutesv.vcf.gz
    """
}
