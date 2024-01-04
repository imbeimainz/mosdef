#' Generate a tidy table with the results of DESeq
#'
#' Generate a tidy table with the results of DESeq
#'
#' @param res_de A \code{\link{DESeqResults}} object
#'
#' @return A "tidy" data.frame with all genes
#' @export
#'
#' @importFrom methods is
#' @importFrom rentrez entrez_summary
#' @importFrom dplyr arrange
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 100, m = 8, betaSD = 1)
#' dds <- DESeq2::DESeq(dds)
#' res <- DESeq2::results(dds)
#' deseqresult2tbl(res)
deseqresult2tbl <- function(res_de) {
  # library("dplyr")
  if (!is(res_de, "DESeqResults")) stop("Not a DESeqResults object.")
  res_de <- as.data.frame(res_de)

  res_de <- cbind(rownames(res_de), res_de)
  names(res_de)[1] <- "id"
  res_de$id <- as.character(res_de$id)

  dplyr::arrange(res_de, .data$padj)
}


geneinfo <- function(gene_id) {
  # the gene id has to be entrez_id

  ## TODO: will need to finish implementation
  entrezinfo <- rentrez::entrez_summary("gene", gene_id)

  return(entrezinfo)
}



