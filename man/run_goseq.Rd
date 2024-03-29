% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_enrichment.R
\name{run_goseq}
\alias{run_goseq}
\title{Extract functional terms enriched in the DE genes, based on goseq}
\usage{
run_goseq(
  res_de = NULL,
  dds = NULL,
  de_genes = NULL,
  bg_genes = NULL,
  top_de = NULL,
  min_counts = 0,
  genome = "hg38",
  id = "ensGene",
  de_type = "up_and_down",
  testCats = c("GO:BP", "GO:MF", "GO:CC"),
  FDR_GO_cutoff = 1,
  mapping = "org.Hs.eg.db",
  add_gene_to_terms = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{res_de}{A DESeqResults object created using \code{DESeq2}}

\item{dds}{A DESeqDataset object created using \code{DESeq2}}

\item{de_genes}{A vector of (differentially expressed) genes}

\item{bg_genes}{A vector of background genes, e.g. all (expressed) genes
in the assays}

\item{top_de}{numeric, how many of the top differentially expressed genes to use for the enrichment analysis.
Attempts to reduce redundancy. Assumes the data is sorted by padj (default in DESeq2).}

\item{min_counts}{numeric, min number of counts a gene needs to have to be included
in the geneset that the de genes are compared to. Default is 0, recommended only for advanced users.}

\item{genome}{A string identifying the genome that genes refer to, as in the
\code{\link{goseq}} function}

\item{id}{A string identifying the gene identifier used by genes, as in the
\code{\link{goseq}} function}

\item{de_type}{One of: 'up', 'down', or 'up_and_down' Which genes to use for GOterm calculations:
upregulated, downregulated or both}

\item{testCats}{A vector specifying which categories to test for over representation amongst DE genes - can be any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG"}

\item{FDR_GO_cutoff}{Numeric value for subsetting the results}

\item{mapping}{Character string, named as the \code{org.XX.eg.db}
package which should be available in Bioconductor}

\item{add_gene_to_terms}{Logical, whether to add a column with all genes annotated
to each GO term}

\item{verbose}{Logical, whether to add messages telling the user which steps were taken}
}
\value{
A table containing the computed GO Terms and related enrichment scores
}
\description{
A wrapper for extracting functional GO terms enriched in a list of (DE) genes,
based on the algorithm and the implementation in the goseq package
}
\details{
Note: the feature length retrieval is based on the \code{\link{goseq}} function,
and requires that the corresponding TxDb packages are installed and available
}
\examples{
library("macrophage")
library("DESeq2")
data(gse, package = "macrophage")

dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage <- DESeq(dds_macrophage)

data(res_de_macrophage, package = "mosdef")
res_de <- res_macrophage_IFNg_vs_naive
mygo <- run_goseq(
  res_de = res_macrophage_IFNg_vs_naive,
  dds = dds_macrophage,
  mapping = "org.Hs.eg.db",
  testCats = "GO:BP",
  add_gene_to_terms = TRUE
)

head(mygo)
}
\seealso{
Other Enrichment functions: 
\code{\link{run_cluPro}()},
\code{\link{run_topGO}()}
}
\concept{Enrichment functions}
