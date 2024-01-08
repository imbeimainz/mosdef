#' Maps numeric values to color values
#'
#' Maps numeric continuous values to values in a color palette
#'
#' @param x A character vector of numeric values (e.g. log2FoldChange values) to
#' be converted to a vector of colors
#' @param pal A vector of characters specifying the definition of colors for the
#' palette, e.g. obtained via \code{\link{brewer.pal}}
#' @param symmetric Logical value, whether to return a palette which is symmetrical
#' with respect to the minimum and maximum values - "respecting" the zero.
#' Defaults to `TRUE`.
#' @param limits A vector containing the limits of the values to be mapped. If
#' not specified, defaults to the range of values in the `x` vector.
#'
#' @return A vector of colors, each corresponding to an element in the original
#' vector
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' a <- 1:9
#' pal <- RColorBrewer::brewer.pal(9, "Set1")
#' map2color(a, pal)
#' plot(a, col = map2color(a, pal), pch = 20, cex = 4)
#'
#' b <- 1:50
#' pal2 <- grDevices::colorRampPalette(
#'   RColorBrewer::brewer.pal(name = "RdYlBu", 11)
#' )(50)
#' plot(b, col = map2color(b, pal2), pch = 20, cex = 3)
map2color <- function(x, pal, symmetric = TRUE, limits = NULL) {
  if (is.null(limits)) {
    limits <- range(x)
  }

  if (symmetric) {
    max_val <- max(limits)
    limits[1] <- -max_val
    limits[2] <- max_val
  }

  pal_ret <- pal[findInterval(x, seq(limits[1],
    limits[2],
    length.out = length(pal) + 1
  ),
  all.inside = TRUE
  )]
  return(pal_ret)
}
