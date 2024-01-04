#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' @param deseqresult A \code{\link{DESeqResults}} object
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#'
#' @return A "tidy" data.frame with only genes marked as differentially expressed
#' @export
#'
#' @importFrom dplyr arrange
#' @importFrom methods is
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 100, m = 8, betaSD = 2)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' deseqresult2DEgenes(res)
deseqresult2DEgenes <- function(deseqresult,
                                FDR = 0.05) {
  # library("dplyr")
  if (!is(deseqresult, "DESeqResults")) stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)

  deseqresult <- cbind(rownames(deseqresult), deseqresult)
  names(deseqresult)[1] <- "id"
  deseqresult$id <- as.character(deseqresult$id)

  # deseqresult$id <- rownames(deseqresult)
  # rownames(deseqresult) <- NULL
  # deseqresult <- dplyr::tbl_df(deseqresult)
  # if("symbol" %in% names(deseqresult))
  #   deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:symbol)
  # else
  #   deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
  tmp <- dplyr::arrange(deseqresult, .data$padj)
  res <- tmp[!(is.na(tmp$padj)) & tmp$padj <= FDR, ]
  res
}



#' Generate a table from the `DESeq2` results
#'
#' Generate a tidy table with the results of `DESeq2`
#'
#' @param res_de A `DESeqResults` object.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to NULL, which would return the full set of results
#' without performing any subsetting based on FDR.
#'
#' @return A tidy `data.frame` with the results from differential expression,
#' sorted by adjusted p-value. If FDR is specified, the table contains only genes
#' with adjusted p-value smaller than the value.
#'
#' @export
#'
#'
#
#' @importFrom methods is
#' 
#' 
#' @examples
#' data(res_de_macrophage, package = "GeneTonic")
#' head(res_macrophage_IFNg_vs_naive)
#' res_df <- deseqresult2df(res_macrophage_IFNg_vs_naive)
#' head(res_df)
deseqresult2df <- function(res_de, FDR = NULL) {
  if (!is(res_de, "DESeqResults")) {
    stop("Not a DESeqResults object.")
  }
  res <- as.data.frame(res_de)
  res <- cbind(rownames(res), res)
  names(res)[1] <- "id"
  res$id <- as.character(res$id)
  res <- res[order(res$padj), ]
  if (!is.null(FDR)) {
    res <- res[!(is.na(res$padj)) & res$padj <= FDR, ]
  }
  res
}
