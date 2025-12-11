#' Generate a table from the `DESeq2` results
#'
#' Generate a tidy table with the results of `DESeq2`
#'
#' @param res_de An object containing the results of the Differential Expression
#' analysis workflow (e.g. `DESeq2`, `edgeR` or `limma`).
#' Currently, this can be a `DESeqResults` object created using the `DESeq2`
#' framework.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to NULL, which would return the full set of
#' results without performing any subsetting based on FDR.
#'
#' @return A tidy `data.frame` with the results from differential expression,
#' sorted by adjusted p-value. If FDR is specified, the table contains only
#' genes with adjusted p-value smaller than the value.
#'
#' @export
#
#' @importFrom methods is
#'
#' @examples
#' library("DESeq2")
#' library("macrophage")
#' data(res_de_macrophage, package = "mosdef")
#' head(res_macrophage_IFNg_vs_naive)
#' res_df <- deresult_to_df(res_macrophage_IFNg_vs_naive)
#' head(res_df)
deresult_to_df <- function(res_de, FDR = NULL) {
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

.deresult_to_df.DESeq <- function() {

}

.deresult_to_df.edgeR <- function() {

}



