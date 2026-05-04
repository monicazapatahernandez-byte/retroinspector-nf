#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Importar módulos
include { MINIMAP2_ALIGN }           from './modules/minimap2'
include { CUTESV }                   from './modules/cutesv'
include { SNIFFLES2 }                from './modules/sniffles2'
include { ASSEMBLY_ALLELES }         from './modules/assembly_alleles'
include { MERGE_INTRAPATIENT }       from './modules/merge_intrapatient'
include { MOSDEPTH_INS; 
          GENOTYPE_INS }             from './modules/genotype_insertions'
include { MERGE_INTERPATIENT }       from './modules/merge_interpatient'
include { REPEATMASKER }             from './modules/repeatmasker'
include { SURVIVOR_INTRAPATIENT;
          MOSDEPTH_DEL;
          GENOTYPE_DEL;
          SURVIVOR_INTERPATIENT;
          GET_DELETIONS }            from './modules/deletions'
include { GET_REFERENCE_REPEATS; GET_REFERENCE_GENOME; GET_REFERENCE_T2T } from './modules/reference'
include { R_PREPARATORY;
          R_GENOTYPING;
          R_ENRICHMENT;
          R_REPORT;
          GENERATE_VCF }            from './modules/r_analysis'

// Mensaje de inicio
log.info """
    R E T R O I N S P E C T O R - N F
    ===================================
    reference   : ${params.reference}
    outdir      : ${params.outdir}
    mode        : ${params.mode}
    callers     : ${params.callers}
    min_reads   : ${params.min_read_support}
    """.stripIndent()

workflow {

    // 1. Canal de muestras desde CSV o directorio
    if (params.input) {
        ch_samples = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, file(row.fastq)) }
    } else if (params.fastq_dir) {
        ch_samples = Channel
            .fromPath("${params.fastq_dir}/*.{fastq,fastq.gz}")
            .map { f -> tuple(f.simpleName, f) }
    } else {
        error "Debes indicar --input (samplesheet CSV) o --fastq_dir"
    }

    // Referencia: usar la proporcionada o descargarla
     if (params.reference) {
        ch_reference = Channel.value(file(params.reference))
    } else if (params.genome == "t2t") {
        GET_REFERENCE_T2T()
        ch_reference = GET_REFERENCE_T2T.out
        log.info "Modo T2T activado -- forzando modo light (sin anotación)"
    } else {
        GET_REFERENCE_GENOME()
        ch_reference = GET_REFERENCE_GENOME.out
    }

    // 2. Alineamiento
    MINIMAP2_ALIGN(ch_samples, ch_reference)
    ch_bams = MINIMAP2_ALIGN.out  // (sample_id, bam, bai)

    // 3. Variant calling
    CUTESV(ch_bams, ch_reference)
    SNIFFLES2(ch_bams, ch_reference)

    // 4. Reconstrucción de secuencias insertadas
    ch_cutesv_input = ch_bams.join(CUTESV.out)
        .map { sid, bam, bai, vcf, csi -> tuple(sid, "cutesv", vcf, csi, bam, bai) }
    ch_sniffles2_input = ch_bams.join(SNIFFLES2.out)
        .map { sid, bam, bai, vcf, csi -> tuple(sid, "sniffles2", vcf, csi, bam, bai) }

    ASSEMBLY_ALLELES(ch_cutesv_input.mix(ch_sniffles2_input))

    // 5. Merge intra-paciente y genotyping de inserciones
    ch_polished = ASSEMBLY_ALLELES.out
        .groupTuple()
        .map { sid, callers, vcfs, csis -> tuple(sid, vcfs, csis) }

    MERGE_INTRAPATIENT(ch_polished)

    ch_for_geno_ins = MERGE_INTRAPATIENT.out.join(ch_bams)
        .map { sid, vcf, bam, bai -> tuple(sid, vcf, bam, bai) }

    MOSDEPTH_INS(ch_for_geno_ins)
    GENOTYPE_INS(MOSDEPTH_INS.out)

    // 6. Merge inter-paciente
    ch_all_vcfs = GENOTYPE_INS.out.map { sid, vcf, csi -> vcf }.collect()
    ch_all_csis = GENOTYPE_INS.out.map { sid, vcf, csi -> csi }.collect()

    MERGE_INTERPATIENT(ch_all_vcfs, ch_all_csis)

    // 7. RepeatMasker
    REPEATMASKER(MERGE_INTERPATIENT.out)

    // 8. Deleciones
    ch_caller_vcfs = CUTESV.out.join(SNIFFLES2.out)
        .map { sid, vcf1, csi1, vcf2, csi2 -> tuple(sid, [vcf1, vcf2], [csi1, csi2]) }

    SURVIVOR_INTRAPATIENT(ch_caller_vcfs)

    ch_for_geno_del = SURVIVOR_INTRAPATIENT.out.join(ch_bams)
        .map { sid, vcf, bam, bai -> tuple(sid, vcf, bam, bai) }

    MOSDEPTH_DEL(ch_for_geno_del)
    GENOTYPE_DEL(MOSDEPTH_DEL.out)

    ch_del_vcfs = GENOTYPE_DEL.out.map { sid, vcf, csi -> vcf }.collect()
    ch_del_csis = GENOTYPE_DEL.out.map { sid, vcf, csi -> csi }.collect()

    SURVIVOR_INTERPATIENT(ch_del_vcfs, ch_del_csis)

    GET_REFERENCE_REPEATS()
    GET_DELETIONS(SURVIVOR_INTERPATIENT.out, GET_REFERENCE_REPEATS.out)

    // 9. Análisis R (solo en modo full)
    ch_sample_ids = ch_samples.map { sid, fastq -> sid }.collect()
    if (params.mode == "full") {
        R_PREPARATORY(REPEATMASKER.out, GET_DELETIONS.out, ch_sample_ids)

        R_GENOTYPING(
            R_PREPARATORY.out[2],  // insertionsTable
            R_PREPARATORY.out[1],   // annotatedInsertionsMin3
            ch_sample_ids
        )

        R_ENRICHMENT(R_GENOTYPING.out[0])  // genes.rds

        GENERATE_VCF(
            R_GENOTYPING.out[1],   // vcf_body
            R_GENOTYPING.out[2],   // vcf_body_lax
            file("${projectDir}/data/header.txt")
        )

        R_REPORT(
            R_PREPARATORY.out[2],  // insertionsTable
            R_PREPARATORY.out[3],  // allIns
            R_PREPARATORY.out[1],  // annotatedInsertionsMin3
            R_PREPARATORY.out[4],  // meDeletionsMin3
            R_ENRICHMENT.out[0],   // egoMF
            R_ENRICHMENT.out[1],   // egoBP
            R_ENRICHMENT.out[2],   // egoCC
            R_ENRICHMENT.out[3],   // do
            R_ENRICHMENT.out[4],   // ncg
            Channel.value([])
        )
    }
}
