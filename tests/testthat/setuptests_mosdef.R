library(airway)
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db")

#Get a dds ad a res_de
data(airway)
airway
dds_airway <- DESeqDataSet(airway, design = ~ cell + dex)
# Example, performing extraction of enriched functional categories in
# detected significantly expressed genes

dds_airway <- DESeq(dds_airway)
res_airway <- results(dds_airway)

# get a vector of de and bg genes
res_subset <- deseqresult2df(res_airway)[1:100, ]
myde <- res_subset$id
myassayed <- rownames(res_airway)