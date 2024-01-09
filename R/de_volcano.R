#' Generates a volcano plot using ggplot2
#' This function generates a base volcanoplot for differentially expressed genes that can then be expanded
#' upon using further ggplot functions.
#'
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param L2FC_cutoff A numeric value that sets the cutoff for the xintercept argument of ggplot
#' @param labeled_genes A numeric value describing the amount of genes to be labeled. This uses the Top(x) highest differntially expressed genes
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#'
#' @return A  \code{ggplot2} volcano plot object that can be extended upon by the user
#' @export
#'
#' @importFrom utils head
#' @importFrom AnnotationDbi mapIds
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline geom_point
#' theme_classic scale_color_manual coord_cartesian scale_x_continuous
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom rlang .data
#'
#'
#'
#' @examples
#' library(ggplot2)
#' library(RColorBrewer) # for a colourful plot
#' library(ggrepel) # for nice annotations
#' library(airway)
#' library(DESeq2)
#' library("org.Hs.eg.db")
#' data(airway)
#' airway
#' dds_airway <- DESeqDataSet(airway, design = ~ cell + dex)
#' # Example, performing extraction of enriched functional categories in
#' # detected significantly expressed genes
#' ## Not run:
#' dds_airway <- DESeq(dds_airway)
#' res_airway <- results(dds_airway)
#' p <- de_volcano(res_airway,
#'   L2FC_cutoff = 1,
#'   labeled_genes = 20,
#'   mapping = "org.Hs.eg.db"
#' )
#' p
de_volcano <- function(res_de,
                       L2FC_cutoff = 1,
                       labeled_genes = 30,
                       mapping = "org.Mm.eg.db") {
  
  if (!is(res_de, "DESeqResults")) {
    stop("The provided `res_de` is not a DESeqResults object, please check your input parameters.")
  }
  
  
  annot_to_map_to <- get(mapping)
  res_de$symbol <- mapIds(annot_to_map_to,
    keys = row.names(res_de),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  df <- deseqresult2df(res_de)

  # finding the highest value in the log2FoldChange column and rounding it up to get a nice symetric plot
  x_limit <- ceiling(max(abs(range(df$log2FoldChange, na.rm = TRUE))))



  df$diffexpressed <- "NO"
  # if log2Foldchange > Cutoff and pvalue < 0.05, set as "UP"
  df$diffexpressed[df$log2FoldChange > L2FC_cutoff & df$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -Cutoff and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df$log2FoldChange < -L2FC_cutoff & df$pvalue < 0.05] <- "DOWN"

  # calculate top 30 degenes based on pvalue
  df$delabel <- ifelse(df$symbol %in% head(df[order(df$pvalue), "symbol"], labeled_genes), df$symbol, NA)




  ggplot(data = df, aes(
    x = .data$log2FoldChange, y = -log10(.data$pvalue),
    colour = .data$diffexpressed, label = .data$delabel
  )) +
    geom_vline(xintercept = c(-L2FC_cutoff, L2FC_cutoff), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    geom_point() +
    theme_classic() +
    scale_color_manual(
      values = c("skyblue", "grey", "tomato"), # to set the colours of our variable
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    coord_cartesian(ylim = c(0, max(-log10(df$pvalue))), xlim = c(-x_limit, x_limit)) +
    scale_x_continuous(breaks = seq(-x_limit, x_limit, 2)) +
    geom_text_repel(max.overlaps = Inf) # To show all labels
}
