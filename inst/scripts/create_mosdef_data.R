library(airway)
library(DESeq2)
library("org.Hs.eg.db")
library("macrophage")
#Get  the base data
data(airway)

# Get a dds object
dds_airway <- DESeqDataSet(airway, design = ~ cell + dex)
dds_airway <- DESeq(dds_airway)

#Get a res_de object
res_airway <- results(dds_airway)

# Get an enrichment object

topgoDE_airway <- topGOtable(
   res_de = res_airway,
   dds = dds_airway,
   ontology = "BP",
   mapping = "org.Hs.eg.db",
   geneID = "symbol",
 )