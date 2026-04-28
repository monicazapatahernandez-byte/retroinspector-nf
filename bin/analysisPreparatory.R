options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(ggplot2)
library(magrittr)
library(qqman)

# === Argumentos de línea de comandos ===
args = commandArgs(trailingOnly = TRUE)
repeatMaskerBedFile    = args[1]
repeatMaskerVcfDelFile = args[2]
minReadSupport         = as.integer(args[3])
allPrefix              = args[4]
samples                = args[5:length(args)]
# =======================================

annotate = function(tableToAnnotate) {
  print("Building annotation")
  library(annotatr)
  library(AnnotationHub)
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  annotationNames = builtin_annotations()[
    stri_detect_fixed(builtin_annotations(), "hg38")
  ]
  annotationNames = annotationNames[annotationNames != "hg38_cpg_inter" &
      annotationNames != "hg38_enhancers_fantom"]
  annotation = build_annotations(genome = "hg38", annotations = annotationNames)
  annotation = annotation[GenomicRanges::mcols(annotation)$type != "hg38_cpg_inter"]
  saveRDS(annotation, "annotation.rds")
  print("Built annotation")
  print("Annotating")
  annotatedInsertions = data.table(data.frame(
      annotatr::annotate_regions(
        regions = GenomicRanges::makeGRangesFromDataFrame(
          tableToAnnotate,
          seqnames.field = "seqnames",
          keep.extra.columns = T),
        annotations = annotation,
        ignore.strand = F,
        quiet = T)))
  print("Annotated")
  annotatedInsertions[, annot.type := stri_replace_first_regex(annot.type, "hg38_(.*)", "$1")]
  text = stri_match_all_regex(annotatedInsertions$annot.type, "([^-]*)_(.*)") %>%
    lapply(., function(x) {x[2:3]})
  annotatedInsertions[, `:=`(annot.class = text %>% lapply(`[`, 1) %>% unlist(),
                 annot.subclass = text %>% lapply(`[`, 2) %>% unlist()
                 )]
  annotatedInsertions[, annot.type := NULL]
  rm(text)
  print("Added class and subclass")
  annotatedInsertions[, c("annot.strand", "annot.seqnames", "annot.width", "annot.id") := NULL]
  print("Finished annotation")
  return(annotatedInsertions)
}

repeatMaskerTable = fread(
    repeatMaskerBedFile,
    col.names = c(
    "seqnames", "start", "end", "name", "score", "strand", "seqId",
    "sw_score",
    "repeat.class",
    "repeat.subclass",
    "repeat.start",
    "repeat.end",
    "repeat.left",
    "repeat.strains",
    "repeat.divergence_percentage",
    "repeat.deletion_percentage",
    "repeat.insertion_percentage",
    "vcf_alt",
    "vcf_info",
    paste0("id", samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))
    )
  )

genoFields = c("GT", "PSV", "LN", "DR", "ST", "QV", "TY", "ID", "RAL", "AAL", "CO", "SC")
genoFieldsOfInt = c("GT", "PSV", "DR")
pos = which(genoFields %in% genoFieldsOfInt)

for (col in repeatMaskerTable %>% colnames() %>%
      .[. %in%
      paste0("id", samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))
      ]) {
  repeatMaskerTable[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
  repeatMaskerTable[, paste0(col, "_DR") := tstrsplit(.SD[[paste0(col, "_DR")]], ",", fixed = T, keep = 2)]
  repeatMaskerTable[,
      paste0(col, "_DR") :=
        fifelse(.SD[[paste0(col, "_DR")]] == ".", 0, as.numeric(.SD[[paste0(col, "_DR")]]) )
      ]
}

repeatMaskerTable[, SUPP_VEC_min3 := apply(.SD, 2, `>=`, y = minReadSupport) %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_DR$", colnames(repeatMaskerTable))]
repeatMaskerTable[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")]

repeatMaskerTable[
  , SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>%
    apply(., 2, as.integer) %>%
    apply(., 1, paste0, collapse = ""),
  .SDcols = patterns("_PSV$")]

bitwiseAndStr = function (a, b) {
  lapply (1:length(a), function (i) {
    a2 = unlist( strsplit(a[[i]], split = "", fixed = T) )
    b2 = unlist( strsplit(b[[i]], split = "", fixed = T) )
    paste0(ifelse(a2 == b2 & a2 == "1", "1", "0"), collapse = "")
  })
}

correctTrusty = unlist( bitwiseAndStr(repeatMaskerTable$SUPP_VEC_min3, repeatMaskerTable$SUPP_VEC_trusty) )
repeatMaskerTable[, SUPP_VEC_trusty := correctTrusty]
rm(correctTrusty)
repeatMaskerTable[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
repeatMaskerTable[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", vcf_info))]

for (col in repeatMaskerTable %>% colnames() %>%
      grep(pattern = "_GT$", x = ., value = T)) {
  repeatMaskerTable[, (col) := get(col) %>%
                                  lapply(function (x) {
                                    x %>%
                                      strsplit(split = "/", fixed = T) %>%
                                      unlist() %>%
                                      as.numeric() %>%
                                      sum()
                                  }) %>% unlist()]
}

repeatMaskerTable[
  , genotype := .SD %>%
    simplify2array() %>%
    t() %>%
    split(., rep(1:ncol(.), each = nrow(.))) %>%
    unname(),
  .SDcols = grep(pattern = "_GT$", colnames(repeatMaskerTable))]
repeatMaskerTable[, SUPP_mask_geno := genotype %>% lapply(`>=`, y = 1)]
repeatMaskerTable[, SUPP_VEC_geno := SUPP_mask_geno %>% lapply(as.numeric) %>% lapply(paste0, collapse = "") %>% unlist()]
repeatMaskerTable[, SUPP_geno := SUPP_mask_geno %>% lapply(sum) %>% unlist()]
repeatMaskerTable[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2 ) ) %>% unlist()]

repeatMaskerTable[, SUPP_mask_min3 := SUPP_VEC_min3 %>%
            lapply(strsplit, split = "") %>%
            lapply(function(x)as.logical(as.numeric(unlist(x))))]
repeatMaskerTable[, SUPP_mask_trusty := SUPP_VEC_trusty %>%
            lapply(strsplit, split = "") %>%
            lapply(function(x)as.logical(as.numeric(unlist(x))))]

repeatMaskerTable[,
  c("repeat.subclass", "repeat.family") := stri_match_all_regex(repeat.subclass, "([^-]+)(?:-([^-]+))?$") %>%
      rapply(., f = `[`, ... = 2:3, how = "r") %>%
      lapply(function (x) {x[is.na(x)] = "-"; x}) %>%
      data.table::transpose()]

repeatMaskerTable[, grep("^id", colnames(repeatMaskerTable)) := NULL]
repeatMaskerTable[, vcf_info := NULL]

chrs = paste0("chr", c(1:22, "X", "Y"))
class1 = c("Retroposon", "SINE", "LINE", "LTR")
class2 = c("DNA")

repeatMaskerTable[, repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN]
repeatMaskerTableMin3 = repeatMaskerTable[SUPP_min3 > 0]
allIns = repeatMaskerTable
repeatMaskerTableMin3 = repeatMaskerTableMin3[repeat.percentage >= 0.85 & repeat.class %in% c(class1, class2) & seqnames %in% chrs]
print("Filtered repeatMaskerTableMin3")
annotatedInsertionsMin3 = annotate(repeatMaskerTableMin3)
print("Removed annotation columns from <insertionsTable>")
saveRDS(annotatedInsertionsMin3, "annotatedInsertionsMin3.rds")
saveRDS(repeatMaskerTableMin3, "insertionsTable.rds")
saveRDS(allIns, "allIns.rds")

# Deletions
meDeletionsMin3 = fread(repeatMaskerVcfDelFile, skip = "#CHROM")
meDeletionsMin3 = meDeletionsMin3[!INFO %like% "ME=[^;]*\\?[^;]*" & !INFO %like% "MEFAM=[^;]*\\?[^;]*"]
meDeletionsMin3 = meDeletionsMin3[`#CHROM` %in% chrs]
meDeletionsMin3[, `:=`(
    survivorId = ID,
    ID = paste0("surv", 1:nrow(meDeletionsMin3)),
    svlen = INFO %>% sub(".*SVLEN=-([0-9]+).*", "\\1", x = .) %>% as.numeric()
  )]

# Columnas de muestra: todo lo que viene después de FORMAT (columna 9)
sample_cols = colnames(meDeletionsMin3)[10:ncol(meDeletionsMin3)]

genoFields = meDeletionsMin3[1, FORMAT] %>% unname() %>% stri_split_fixed(pattern = ":") %>% unlist()
genoFieldsOfInt = c("GT", "PSV")
pos = which(genoFields %in% genoFieldsOfInt)
for (col in sample_cols) {
  meDeletionsMin3[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
}
meDeletionsMin3[, SUPP_VEC_min3 := apply(.SD, 2, `!=`, y = "NaN") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_PSV$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")]
meDeletionsMin3[, SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = patterns("_PSV$")]
meDeletionsMin3[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
meDeletionsMin3[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", INFO)) %>% abs()]

for (col in meDeletionsMin3 %>% colnames() %>%
     grep(pattern = "_GT$", x = ., value = T)) {
  meDeletionsMin3[, (col) := get(col) %>% ifelse(. == "./.", "0/0", .) %>%
                    lapply(function (x) {
                      x %>%
                        strsplit(split = "/", fixed = T) %>%
                        unlist() %>%
                        as.numeric() %>%
                        sum()
                    }) %>% unlist()]
}

meDeletionsMin3[
  , genotype := .SD %>%
    simplify2array() %>%
    t() %>%
    split(., rep(1:ncol(.), each = nrow(.))) %>%
    unname(),
  .SDcols = grep(pattern = "_GT$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_VEC_geno := genotype %>% lapply(`>=`, y = 1)]
meDeletionsMin3[, SUPP_geno := SUPP_VEC_geno %>% lapply(sum) %>% unlist()]
meDeletionsMin3[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2) ) %>% unlist()]
meDeletionsMin3[, c("FORMAT", sample_cols) := NULL]
meDeletionsMin3 = meDeletionsMin3[
  , .(
    seqnames = `#CHROM`,
    start = POS,
    end = sub(pattern = ".*;END=([^;]+).*", replacement = "\\1", x = INFO, perl = T) %>% as.numeric(),
    id = ID,
    SUPP_min3, SUPP_VEC_min3, SUPP_trusty, SUPP_VEC_trusty,
    genotype, maf, SUPP_geno, SUPP_VEC_geno, svlen,
    name = sub(pattern = ".*ME=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.class = sub(pattern = ".*MEFAM=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.coords = sub(pattern = ".*MECOORDS=([^;]+).*", replacement = "\\1", x = INFO, perl = T)
  )
][, c("repeat.class", "repeat.subclass") := repeat.class %>%
    stri_match_all_regex("([^/]+)(?:/(.+))?$") %>%
    rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()][
      , repeat.subclass := fifelse(is.na(repeat.subclass), "-", repeat.subclass)
    ][
      , c("repeat.subclass", "repeat.family") :=
        stri_match_all_regex(repeat.subclass, "([^-]+|-)(?:-(.+))?$") %>%
        rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()
    ][, repeat.family := fifelse(is.na(repeat.family), "-", repeat.family)]

saveRDS(meDeletionsMin3, "meDeletionsMin3.rds")
