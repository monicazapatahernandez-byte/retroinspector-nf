library(data.table)
library(magrittr)
library(ggplot2)
library(patchwork)
library(stringi)
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(DOSE)

# === Argumentos de línea de comandos ===
args = commandArgs(trailingOnly = TRUE)
genesFile = args[1]
p         = as.numeric(args[2])
# =======================================

genes     = readRDS(genesFile)
geneNames = genes$annot.symbol %>% unique()
geneDt    = bitr(geneNames, fromType = "SYMBOL",
                 toType = c("REFSEQ", "ENSEMBL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db) %>% data.table()

# GO enrichment
for (x in c("MF", "BP", "CC")) {
  name = paste0("ego", x)
  assign(x = name,
         enrichGO(
           gene          = geneDt$ENSEMBL,
           OrgDb         = org.Hs.eg.db,
           keyType       = 'ENSEMBL',
           ont           = x,
           pAdjustMethod = "BH",
           pvalueCutoff  = p)
  )
  object = get(name)
  if (!is.null(object)) {
    object = object %>% setReadable(OrgDb = org.Hs.eg.db)
    assign(name, object)
  }
  saveRDS(get(name), paste0(name, ".rds"))
}

# DO
do = enrichDO(gene          = geneDt$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = p,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              readable      = T)
saveRDS(do, "do.rds")

# NCG
ncg = enrichNCG(geneDt$ENTREZID %>% unique(),
                pvalueCutoff = p,
                readable     = T)
saveRDS(ncg, "ncg.rds")
