process SNIFFLES2 {
    tag "$sample_id"
    conda "${projectDir}/env_sniffles2.yaml"

    publishDir "${params.outdir}/variants/sniffles2", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.sniffles2.vcf.gz"), path("${sample_id}.sniffles2.vcf.gz.csi")

    script:
    def fixScript = "${projectDir}/bin/fixVCF.py"
    """
    sniffles \
        --threads ${task.cpus} \
        --output-rnames \
        --input ${bam} \
        --vcf ${sample_id}.sniffles2.temp.vcf

    cat ${sample_id}.sniffles2.temp.vcf | \
    gawk -v 'OFS=\t' \
    '{if (substr(\$0,1,1)=="#") {print} \
    else { match(\$8,/SVTYPE=([^;]+)/,a); \
      if (a[1]=="INS"||a[1]=="BND") {print} \
      else {\$5="<"a[1]">";\$4="N";print} }}' | \
    gawk -v 'OFS=\t' \
    '{if (substr(\$0,1,1)=="#") {print} \
    else { match(\$8,/SUPPORT=([0-9]+)/,a); \
      if (a[1]+0 >= ${params.min_read_support}) {print} }}' | \
    gawk -v 'OFS=\t' \
    '{if (substr(\$0,1,1)=="#") {sub(/SUPPORT,/,"RE,",\$0);print} \
    else {sub(/SUPPORT=/,"RE=",\$8);print}}' | \
    python3 ${fixScript} | \
    bcftools sort -O z -o ${sample_id}.sniffles2.vcf.gz

    bcftools index ${sample_id}.sniffles2.vcf.gz
    """
}
