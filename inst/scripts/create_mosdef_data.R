library("GeneTonic")
library("macrophage")
data(gse)

# dds object -------------------------------------------------------------------
library("DESeq2")
dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
# no need to save this one, can be readily generated

# res object -------------------------------------------------------------------
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
library("org.Hs.eg.db")
dds_macrophage <- DESeq(dds_macrophage)

res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.05)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL
library("AnnotationDbi")

save(res_macrophage_IFNg_vs_naive, file = "data/res_de_macrophage.RData", compress = "xz")

# res_enrich object topGO ------------------------------------------------------
library("topGO")
topgoDE_macrophage_IFNg_vs_naive <-topGOtable(
  res_de = res_macrophage_IFNg_vs_naive,
  dds = dds_macrophage,
  ontology = "BP",
  mapping = "org.Hs.eg.db",
  geneID = "symbol",
)

save(topgoDE_macrophage_IFNg_vs_naive, file = "data/res_enrich_macrophage_topGO.RData", compress = "xz")

# res_enrich object goseq ------------------------------------------------------
library("goseq")
goseqDE_macrophage_IFNg_vs_naive <-goseqTable(
  res_de = res_macrophage_IFNg_vs_naive,
  dds = dds_macrophage,
  mapping = "org.Hs.eg.db",
  testCats = "GO:BP",
  add_gene_to_terms = TRUE
)
save(goseqDE_macrophage_IFNg_vs_naive, file = "data/res_enrich_macrophage_goseq.RData", compress = "xz")

# res_enrich object clusterProfiler --------------------------------------------
library(clusterProfiler)
cluProDE_macrophage_IFNg_vs_naive <- cluproTable(
  res_de = res_macrophage_IFNg_vs_naive,
  dds = dds_macrophage,
  mapping = "org.Hs.eg.db"
)
save(topgoDE_macrophage_IFNg_vs_naive, file = "data/res_enrich_macrophage_cluPro.RData", compress = "xz")

