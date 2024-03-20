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


#' Style DT color bars
#'
#' Style DT color bars for values that diverge from 0.
#'
#' @details This function draws background color bars behind table cells in a column,
#' width the width of bars being proportional to the column values *and* the color
#' dependent on the sign of the value.
#'
#' A typical usage is for values such as `log2FoldChange` for tables resulting from
#' differential expression analysis.
#' Still, the functionality of this can be quickly generalized to other cases -
#' see in the examples.
#'
#' The code of this function is heavily inspired from styleColorBar, and borrows
#' at full hands from an excellent post on StackOverflow -
#' https://stackoverflow.com/questions/33521828/stylecolorbar-center-and-shift-left-right-dependent-on-sign/33524422#33524422
#'
#' @param data The numeric vector whose range will be used for scaling the table
#' data from 0-100 before being represented as color bars. A vector of length 2
#' is acceptable here for specifying a range possibly wider or narrower than the
#' range of the table data itself.
#' @param color_pos The color of the bars for the positive values
#' @param color_neg The color of the bars for the negative values
#'
#' @return This function generates JavaScript and CSS code from the values
#' specified in R, to be used in DT tables formatting.
#'
#' @export
#'
#' @importFrom DT JS
#'
#' @examples
#'
#' # With a very simple data frame
#'
#' simplest_df <- data.frame(
#'   a = c(rep("a", 9)),
#'   value = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
#' )
#'
#' library("DT")
#' DT::datatable(simplest_df) |>
#'   formatStyle(
#'     "value",
#'     background = styleColorBar_divergent(
#'       simplest_df$value,
#'       scales::alpha("forestgreen", 0.4),
#'       scales::alpha("gold", 0.4)
#'     ),
#'     backgroundSize = "100% 90%",
#'     backgroundRepeat = "no-repeat",
#'     backgroundPosition = "center"
#'   )
styleColorBar_divergent <- function(data,
                                    color_pos,
                                    color_neg) {
  max_val <- max(abs(data))
  JS(
    sprintf(
      "isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
      max_val, color_pos, max_val, color_pos, color_neg, color_neg, max_val, max_val
    )
  )
}


#' DE table painter
#'
#' Beautifying the aspect and looks of a DE results table
#'
#' Feeding on the classical results of DE workflows, this function formats and
#' tries to prettify the representation of the key values in it.
#'
#' @param res_de A `DESeqResults` object created using `DESeq2`, or a data frame
#' obtained from such an object through [deseqresult2df()]
#' @param rounding_digits Numeric value, specifying the number of digits to round
#' the numeric values of the DE table (except the p-values)
#' @param signif_digits Numeric value, specifying the number of significant digits
#' to display for the p-values in the DE table
#' @param up_DE_color Character string, specifying the color to use for coloring
#' the bar of upregulated genes.
#' @param down_DE_color Character string, specifying the color to use for coloring
#' the bar of downregulated genes.
#' @param logfc_column Character string, defining the name of the column in which
#' to find the log2 fold change.
#' @param basemean_column Character string, defining the name of the column in which
#' to find the average expression value.
#' @param lfcse_column Character string, defining the name of the column in which
#' to find the standard error of the log2 fold change.
#' @param stat_column Character string, defining the name of the column in which
#' to find the values of the test statistic.
#' @param pvalue_column Character string, defining the name of the column in which
#' to find the unadjusted p-values.
#' @param padj_column Character string, defining the name of the column in which
#' to find the adjusted p-values.
#'
#' @return A `datatable` object, ready to be rendered as a widget inside an
#' analysis Rmarkdown report.
#'
#' @export
#'
#' @importFrom scales alpha
#' @importFrom DT datatable formatRound formatSignif formatStyle
#'
#' @examples
#' data(res_de_macrophage, package = "mosdef")
#' de_table_painter(res_macrophage_IFNg_vs_naive,
#'                  rounding_digits = 3,
#'                  signif_digits = 5)
#'
#' ## It is also possible to pass the "buttonified" table,
#' res_df_small <- deseqresult2df(res_macrophage_IFNg_vs_naive)[1:100, ]
#'
#' buttonified_df <- buttonifier(res_df_small,
#'                               new_cols = c("NCBI", "HPA"),
#'                               ens_col = "id",
#'                               ens_species = "Homo_sapiens",
#'                               output_format = "DF"
#' )
#'
#' de_table_painter(buttonified_df,
#'                  rounding_digits = 3,
#'                  signif_digits = 5)
de_table_painter <- function(res_de,
                             rounding_digits = NULL,
                             signif_digits = NULL,
                             up_DE_color = "darkred",
                             down_DE_color = "navyblue",
                             logfc_column = "log2FoldChange",
                             basemean_column = "baseMean",
                             lfcse_column = "lfcSE",
                             stat_column = "stat",
                             pvalue_column = "pvalue",
                             padj_column = "padj"
                             ) {
  ## Checks on the input parameters

  res_de <- res_de[order(res_de$padj), ]
  my_dt <- DT::datatable(
    as.data.frame(res_de),
    escape = FALSE,
    options = list(
      scrollX = TRUE,
      scrollY = "400px",
      pageLength = 25,
      columnDefs = list(
        list(className = "dt-center", targets = "_all")
      )
    )
  )

  if (!is.null(rounding_digits)) {
    my_dt <- formatRound(table = my_dt,
                         columns = c(logfc_column, basemean_column, lfcse_column, stat_column),
                         digits = rounding_digits)
  }

  if (!is.null(signif_digits)) {
    my_dt <- formatSignif(table = my_dt,
                          columns = c(pvalue_column, padj_column),
                          digits = signif_digits)
  }

  my_dt <- formatStyle(
    table = my_dt,
    columns = logfc_column,
    background = styleColorBar_divergent(
      as.data.frame(res_de)[[logfc_column]],
      scales::alpha(down_DE_color, 0.4),
      scales::alpha(up_DE_color, 0.4)
    ),
    backgroundSize = "100% 90%",
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center"
  )

  return(my_dt)
}
