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
insertionsTableFile      = args[1]
annotatedInsertionsFile  = args[2]
allPrefix                = args[3]
samples                  = allPrefix
# =======================================

chrs   = paste0("chr", c(1:22, "X", "Y"))
class1 = c("Retroposon", "SINE", "LINE", "LTR")
class2 = c("DNA")

insertionsTable       = readRDS(insertionsTableFile)
annotatedInsertionsMin3 = readRDS(annotatedInsertionsFile)

# Genes
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

print("Processing deletions")

# Write VCF output
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

saveRDS(genes, "genes.rds")
