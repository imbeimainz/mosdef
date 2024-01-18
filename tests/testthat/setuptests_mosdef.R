library(airway)
library(DESeq2)
library(topGO)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("macrophage")

#Get a dds ad a res_de
data(airway)
airway
dds_airway_nodeseq <- DESeqDataSet(airway, design = ~ cell + dex)
# Example, performing extraction of enriched functional categories in
# detected significantly expressed genes

dds_airway <- DESeq(dds_airway_nodeseq)
res_airway <- results(dds_airway)
res_airway_nosymbols <- results(dds_airway)

res_airway$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys = row.names(res_airway),
                                       column = "SYMBOL",
                                       keytype = "ENSEMBL",
                                       multiVals = "first"
)
res_airway$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                           keys = row.names(res_airway),
                                           column = "SYMBOL",
                                           keytype = "ENSEMBL",
                                           multiVals = "first"
)
airway_df <- deseqresult2df(res_airway)

# get a vector of de and bg genes
res_subset <- deseqresult2df(res_airway)[1:500,]
myde <- res_subset$id
myassayed <- rownames(res_airway)
annotationobject <- deseqresult2df(res_airway)
annotationobject <- annotationobject["SYMBOL"]
#Macrophage


# dds --------------------------------------------------------------------------
data(gse)
dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

# annotation -------------------------------------------------------------------
anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db,
                     keys = rownames(dds_macrophage),
                     column = "SYMBOL",
                     keytype = "ENSEMBL"
  ),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)
# alternatively, one could use the wrapper in ...
# anno_df <- pcaExplorer::get_annotation_orgdb(dds_macrophage, "org.Hs.eg.db", "ENSEMBL")

# res_de -----------------------------------------------------------------------
## using counts and average transcript lengths from tximeta
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_unnormalized <- dds_macrophage

dds_macrophage <- DESeq(dds_macrophage)
vst_macrophage <- vst(dds_macrophage)
res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.05
)
summary(res_macrophage_IFNg_vs_naive)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL


goseqde_airway <- goseqTable(
  res_de = res_airway,
  dds = dds_airway,
  mapping = "org.Hs.eg.db"
)