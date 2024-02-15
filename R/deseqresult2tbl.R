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
#' @importFrom rlang .data
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' data(res_airway, package = "mosdef")
#' deseqresult2tbl(res_airway)
deseqresult2tbl <- function(res_de) {
  # library("dplyr")
  if (!is(res_de, "DESeqResults")) stop("Not a DESeqResults object.")
  res_de <- as.data.frame(res_de)

  res_de <- cbind(rownames(res_de), res_de)
  names(res_de)[1] <- "id"
  res_de$id <- as.character(res_de$id)

  dplyr::arrange(res_de, .data$padj)
}

