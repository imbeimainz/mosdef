#' Pairwise scatter plot matrix and correlation plot of counts
#'
#' @param df A data frame, containing the (raw/normalized/transformed) counts
#' @param log Logical, whether to convert the input values to log2 (with addition
#' of a pseudocount). Defaults to TRUE.
#' @param method Character string, one of `pearson` (default), `kendall`, or
#' `spearman` as in `cor`
#' @param use_subset Logical value. If TRUE, only 1000 values per sample will be used
#' to speed up the plotting operations.
#'
#' @return A plot with pairwise scatter plots and correlation coefficients
#'
#' @export
#'
#' @importFrom grDevices colorRamp rgb
#' @importFrom graphics pairs par strwidth text
#' @importFrom stats cor
#'
#' @examples
#'
#' library("macrophage")
#' library("DESeq2")
#' data(gse, package = "macrophage")
#' ## dds object
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' ## Using just a subset for the example
#' pair_corr(counts(dds_macrophage, normalized = TRUE)[1:100, 1:8])
pair_corr <- function(df, log = TRUE, method = "pearson", use_subset = TRUE) {
  if (log) {
    df <- log2(1 + df)
  }

  if (use_subset) {
    set.seed(42)
    df <- df[sample(1:nrow(df),
                    min(nrow(df), 1000)), ]
  }

  ## Get min and max count values for axis range.
  rangeMin <- min(df)
  rangeMax <- max(df)

  colorFunction <- colorRamp(c("black", "red"))
  ## colorFunction() expects values from 0 to 1.
  zMatrix <- colorFunction(seq(0, 1, by = .01))
  # zColors goes from 1 to 100.
  zColors <- sort(rgb(zMatrix[, 1], zMatrix[, 2], zMatrix[, 3], maxColorValue = 255))
  labelSize <- 1
  title <- "Pairwise Correlations"
  ## Modified from R pairs() documentation
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method = method))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")

    # color text based on r value and change size of text also based on r value (larger text for larger r value).
    cex.cor <- labelSize / strwidth(txt)
    # color text based on r value (red is r=1).
    text(0.5, 0.5, txt, cex = cex.cor * r * 0.7, col = zColors[r * 100])
  }

  pairs(df, pch = 20, col = alpha("black", 0.4),
        cex.labels = labelSize,
        main = title,
        upper.panel = panel.cor,
        xlim = c(rangeMin, rangeMax),
        ylim = c(rangeMin, rangeMax))
}
