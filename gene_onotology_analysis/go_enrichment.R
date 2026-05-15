library(clusterProfiler)
library(org.Hs.eg.db)

genes <- read.table(
  "genes_with_motif.txt",
  stringsAsFactors=FALSE
)[,1]

genes <- sub("::.*", "", genes)

entrez <- bitr(
  genes,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)

ego <- enrichGO(
  gene          = entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

pdf("go_enrichment_results.pdf", width=12, height=7)

dotplot(ego, showCategory=15)

dev.off()
