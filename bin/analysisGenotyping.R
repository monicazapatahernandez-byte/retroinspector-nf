options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
print("Loading libraries")
library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(magrittr)
library(Rcpp)
library(parallel)
print("Loaded libraries")
# === Argumentos de línea de comandos ===
args = commandArgs(trailingOnly = TRUE)
insertionsTableFile = args[1]
annotatedInsertionsFile = args[2]
allPrefix = args[3]
meDeletionsFile = args[4]
samples = args[5:length(args)]
# =======================================
chrs   = paste0("chr", c(1:22, "X", "Y"))
class1 = c("Retroposon", "SINE", "LINE", "LTR")
class2 = c("DNA")
insertionsTable       = readRDS(insertionsTableFile)
annotatedInsertionsMin3 = readRDS(annotatedInsertionsFile)
print("Processing gene info")
genes = annotatedInsertionsMin3[
  SUPP_trusty > 0 & SUPP_geno > 0 &
  annot.subclass %in% c("cds", "introns", "promoters") &
  !is.na(annot.symbol),
  .SD[1],
  by = c("seqId", "annot.symbol")][order(annot.symbol)]
genes = genes[, .(
  uniques     = seqId %>% unique() %>% length(),
  effectLax   = sum(SUPP_min3),
  effectStrict = sum(SUPP_trusty)
), by = annot.symbol][order(-effectStrict)]
print("Finished gene processing")
roi = annotatedInsertionsMin3[
  repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") &
  SUPP_trusty > 0 &
  repeat.percentage >= 0.85]
ioi2 = roi[, .SD[1], by = seqId][,
  .(seqId, seqnames, start, end,
    name, repeat.class, repeat.subclass,
    SVLEN, SUPP_min3, SUPP_trusty,
    alleles = genotype %>% lapply(sum) %>% unlist()
  )]
setkey(ioi2, seqnames, start, end)
print("Processing insertions")
print("Writing VCF")
vcf = ioi2[insertionsTable[, c("seqId", "SUPP_geno", "maf", "vcf_alt")], on = "seqId", nomatch = NULL]
vcf = vcf[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS      = start,
  ID       = seqId,
  REF      = "N",
  ALT      = vcf_alt,
  QUAL     = ".",
  FILTER   = ".",
  INFO     = paste(
    "SVTYPE=INS",
    paste0("TE_FAMILY=",    repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=",      name),
    paste0("SVLEN=",        SVLEN),
    paste0("END=",          start),
    paste0("SUPP_strict=",  SUPP_trusty),
    paste0("SUPP_lax=",     SUPP_min3),
    paste0("SUPP_geno=",    SUPP_geno),
    paste0("MAF=",          maf),
    paste0("SAMPLE_N=",     length(samples)),
    paste0("ALLELE_N=",     alleles),
    sep = ";")
)]
fwrite(vcf, paste0(allPrefix, ".me.insertions.txt"), sep = "\t")

# Exportar input para hallmarks.py (inserciones estrictas)
fwrite(
    vcf[, .(seqId = ID, seqnames = `#CHROM`, start = POS, vcf_alt = ALT)],
    paste0(allPrefix, ".hallmarks_input.tsv"),
    sep = "\t"
)
ioi3 = annotatedInsertionsMin3[
  repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") &
  repeat.percentage >= 0.85][, .SD[1], by = seqId][,
  .(seqId, seqnames, start, end,
    name, repeat.class, repeat.subclass,
    SVLEN, SUPP_min3, SUPP_trusty,
    alleles = genotype %>% lapply(sum) %>% unlist()
  )]
setkey(ioi3, seqnames, start, end)
vcf2 = ioi3[insertionsTable[, c("seqId", "SUPP_geno", "maf", "vcf_alt")], on = "seqId", nomatch = NULL]
vcf2 = vcf2[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS      = start,
  ID       = seqId,
  REF      = "N",
  ALT      = vcf_alt,
  QUAL     = ".",
  FILTER   = ".",
  INFO     = paste(
    "SVTYPE=INS",
    paste0("TE_FAMILY=",    repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=",      name),
    paste0("SVLEN=",        SVLEN),
    paste0("END=",          start),
    paste0("SUPP_strict=",  SUPP_trusty),
    paste0("SUPP_lax=",     SUPP_min3),
    paste0("SUPP_geno=",    SUPP_geno),
    paste0("MAF=",          maf),
    paste0("SAMPLE_N=",     length(samples)),
    paste0("ALLELE_N=",     alleles),
    sep = ";")
)]
fwrite(vcf2, paste0(allPrefix, ".me.insertions.lax.txt"), sep = "\t")
print("Processing deletions")
del = readRDS(meDeletionsFile)
del_roi = del[
  repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") &
  SUPP_trusty > 0]
del_ioi2 = del_roi[, .SD[1], by = id][,
  .(id, seqnames, start, end,
    name, repeat.class, repeat.subclass,
    svlen, SUPP_min3, SUPP_trusty,
    alleles = genotype %>% lapply(sum) %>% unlist()
  )]
setkey(del_ioi2, seqnames, start, end)
vcf_del = del_ioi2[del[, c("id", "SUPP_geno", "maf")], on = "id", nomatch = NULL]
vcf_del = vcf_del[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS      = start,
  ID       = id,
  REF      = "N",
  ALT      = "<DEL>",
  QUAL     = ".",
  FILTER   = ".",
  INFO     = paste(
    "SVTYPE=DEL",
    paste0("TE_FAMILY=",    repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=",      name),
    paste0("SVLEN=-",       svlen),
    paste0("END=",          end),
    paste0("SUPP_strict=",  SUPP_trusty),
    paste0("SUPP_lax=",     SUPP_min3),
    paste0("SUPP_geno=",    SUPP_geno),
    paste0("MAF=",          maf),
    paste0("SAMPLE_N=",     length(samples)),
    paste0("ALLELE_N=",     alleles),
    sep = ";")
)]
fwrite(vcf_del, paste0(allPrefix, ".me.deletions.txt"), sep = "\t")
del_ioi3 = del[
  repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR")][, .SD[1], by = id][,
  .(id, seqnames, start, end,
    name, repeat.class, repeat.subclass,
    svlen, SUPP_min3, SUPP_trusty,
    alleles = genotype %>% lapply(sum) %>% unlist()
  )]
setkey(del_ioi3, seqnames, start, end)
vcf_del2 = del_ioi3[del[, c("id", "SUPP_geno", "maf")], on = "id", nomatch = NULL]
vcf_del2 = vcf_del2[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS      = start,
  ID       = id,
  REF      = "N",
  ALT      = "<DEL>",
  QUAL     = ".",
  FILTER   = ".",
  INFO     = paste(
    "SVTYPE=DEL",
    paste0("TE_FAMILY=",    repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=",      name),
    paste0("SVLEN=-",       svlen),
    paste0("END=",          end),
    paste0("SUPP_strict=",  SUPP_trusty),
    paste0("SUPP_lax=",     SUPP_min3),
    paste0("SUPP_geno=",    SUPP_geno),
    paste0("MAF=",          maf),
    paste0("SAMPLE_N=",     length(samples)),
    paste0("ALLELE_N=",     alleles),
    sep = ";")
)]
fwrite(vcf_del2, paste0(allPrefix, ".me.deletions.lax.txt"), sep = "\t")
genes_del = del[
  SUPP_trusty > 0 & SUPP_geno > 0 &
  !is.na(name),
  .SD[1],
  by = c("id", "name")][order(name)]
genes_del = genes_del[, .(
  uniques     = id %>% unique() %>% length(),
  effectLax   = sum(SUPP_min3),
  effectStrict = sum(SUPP_trusty)
), by = name][order(-effectStrict)]
saveRDS(genes_del, paste0(allPrefix, ".genes_deletions.rds"))
saveRDS(genes, "genes.rds")
