library(DESeq2)
library(topGO)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("macrophage")
library("scales")
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")

data(gse)
dds_macrophage_nodeseq <- DESeqDataSet(gse, design = ~line + condition)
# Example, performing extraction of enriched functional categories in
# detected significantly expressed genes

dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]

# res_de

data(res_de_macrophage, package = "mosdef")

res_macrophage_IFNg_vs_naive$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                           keys = row.names(res_macrophage_IFNg_vs_naive),
                                           column = "SYMBOL",
                                           keytype = "ENSEMBL",
                                           multiVals = "first"
)
res_macrophage_IFNg_vs_naive$symbol <- res_macrophage_IFNg_vs_naive$SYMBOL
macrophage_df <- deseqresult2df(res_macrophage_IFNg_vs_naive)

# get a vector of de and bg genes
res_subset <- deseqresult2df(res_macrophage_IFNg_vs_naive)[1:500, ]
myde <- res_subset$id
myassayed <- rownames(res_macrophage_IFNg_vs_naive)
annotationobject <- deseqresult2df(res_macrophage_IFNg_vs_naive)
annotationobject <- annotationobject["SYMBOL"]


# airway 
library(airway)
#Get  the base data
data(airway)

# Get a dds object and a res object
dds_macrophage <- DESeqDataSet(airway, design = ~ cell + dex)
dds_macrophage <- DESeq(dds_macrophage)

res_macrophage_IFNg_vs_naive_nosymbols <- results(dds_macrophage)

# res_enrich
data(res_enrich_macrophage_topGO, package = "mosdef")



