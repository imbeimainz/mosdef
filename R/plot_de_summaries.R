#' MA-plot from base means and log fold changes
#'
#' MA-plot from base means and log fold changes, in the ggplot2 framework, with
#' additional support to annotate genes if provided.
#'
#' The genes of interest are to be provided as gene symbols if a \code{symbol}
#' column is provided in \code{res_de}, or else by using  the identifiers specified
#' in the row names
#'
#' @param res_de A \code{\link{DESeqResults}} object
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#' @param point_alpha Alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param sig_color Color to use to mark differentially expressed genes. Defaults to red
#' @param annotation_obj A \code{data.frame} object, with row.names as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. Optional
#' @param draw_y0 Logical, whether to draw the horizontal line at y=0. Defaults to
#' TRUE.
#' @param hlines The y coordinate (in absolute value) where to draw horizontal lines,
#' optional
#' @param title A title for the plot, optional
#' @param xlab X axis label, defaults to "mean of normalized counts - log10 scale"
#' @param ylim Vector of two numeric values, Y axis limits to restrict the view
#' @param add_rug Logical, whether to add rug plots in the margins
#' @param intgenes Vector of genes of interest. Gene symbols if a \code{symbol}
#' column is provided in \code{res_de}, or else the identifiers specified in the
#' row names
#' @param intgenes_color The color to use to mark the genes on the main plot.
#' @param labels_intgenes Logical, whether to add the gene identifiers/names close
#' to the marked plots
#' @param labels_repel Logical, whether to use \code{geom_text_repel} for placing the
#' labels on the features to mark
#'
#' @return An object created by \code{ggplot}
#' @export
#'
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_point geom_text geom_rug
#'  scale_colour_manual coord_cartesian xlab ylab ggtitle theme_bw
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#'
#'
#' @examples
#' data(res_de_macrophage, package = "mosdef")
#' plot_ma(res_macrophage_IFNg_vs_naive, FDR = 0.05, hlines = 1)
#'
#' plot_ma(res_macrophage_IFNg_vs_naive,
#'   FDR = 0.1,
#'   intgenes = c(
#'     "ENSG00000103196", # CRISPLD2
#'     "ENSG00000120129", # DUSP1
#'     "ENSG00000163884", # KLF15
#'     "ENSG00000179094" # PER1
#'   )
#' )
plot_ma <- function(res_de,
                    FDR = 0.05,
                    point_alpha = 0.2,
                    sig_color = "red",
                    annotation_obj = NULL, # TODO: add a check, if not available skip this part
                    draw_y0 = TRUE,
                    hlines = NULL,
                    title = NULL,
                    xlab = "mean of normalized counts - log10 scale",
                    ylim = NULL,
                    add_rug = TRUE,
                    intgenes = NULL,
                    intgenes_color = "steelblue",
                    labels_intgenes = TRUE,
                    labels_repel = TRUE) {
  if (!is(res_de, "DESeqResults")) {
    stop("The provided `res_de` is not a DESeqResults object, please check your input parameters.")
  }



  ma_df <- data.frame(
    mean = res_de$baseMean,
    lfc = res_de$log2FoldChange,
    padj = res_de$padj,
    isDE = ifelse(is.na(res_de$padj), FALSE, res_de$padj < FDR),
    ID = rownames(res_de)
  )

  ma_df <- ma_df[ma_df$mean > 0, ]

  if (!is.null(annotation_obj)) {
    ma_df$genename <- annotation_obj$gene_name[match(ma_df$ID, rownames(annotation_obj))]
  }

  ma_df$logmean <- log10(ma_df$mean) # TO ALLOW FOR BRUSHING!!
  # ma_df$DE <- ifelse(ma_df$isDE,"yes","no")
  ma_df$DE <- ifelse(ma_df$isDE, "red", "black")

  p <- ggplot(ma_df, aes(x = .data$logmean, y = .data$lfc, colour = .data$DE))

  if (!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) +
      geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
  }

  if (draw_y0) {
    p <- p + geom_hline(aes(yintercept = 0), col = "red", alpha = 0.4)
  }

  p <- p + xlab(xlab) + ylab("log fold change")

  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(
    name = paste0("FDR = ", FDR),
    values = c("black", sig_color),
    labels = c("nonDE", "DE")
  )

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  if (!is.null(intgenes)) {
    # now here for the symbol
    res_df <- as.data.frame(res_de)
    res_df$logmean <- log10(res_df$baseMean)

    if ("symbol" %in% colnames(res_df)) {
      # use the gene names
      df_intgenes <- res_df[res_df$symbol %in% intgenes, ]
      df_intgenes$myids <- df_intgenes$symbol
    } else {
      # use whatever is there as id
      df_intgenes <- res_df[rownames(res_df) %in% intgenes, ]
      df_intgenes$myids <- rownames(df_intgenes)
    }

    # df_intgenes <- res_df[res_df$symbol %in% intgenes,]
    p <- p + geom_point(data = df_intgenes, aes(.data$logmean, .data$log2FoldChange), color = intgenes_color, size = 4)

    if (labels_intgenes) {
      if (labels_repel) {
        p <- p + geom_text_repel(
          data = df_intgenes, aes(.data$logmean, .data$log2FoldChange, label = .data$myids),
          color = intgenes_color, size = 5
        )
      } else {
        p <- p + geom_text(
          data = df_intgenes, aes(.data$logmean, .data$log2FoldChange, label = .data$myids),
          color = intgenes_color, size = 5, hjust = 0.25, vjust = -0.75
        )
      }
    }
  }

  if (add_rug) {
    p <- p + geom_rug(alpha = 0.3)
  }

  p <- p + theme_bw()
  p
}


#' Generates a volcano plot using ggplot2
#' This function generates a base volcanoplot for differentially expressed genes that can then be expanded
#' upon using further ggplot functions.
#'
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param L2FC_cutoff A numeric value that sets the cutoff for the xintercept argument of ggplot
#' @param labeled_genes A numeric value describing the amount of genes to be labeled. This uses the Top(x) highest
#' differentially expressed genes
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#'
#' @return A  \code{ggplot2} volcano plot object that can be extended upon by the user
#' @export
#'
#' @importFrom methods is
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
#' library(DESeq2)
#' library("org.Hs.eg.db")
#' data(res_de_macrophage, package = "mosdef")
#' p <- de_volcano(res_macrophage_IFNg_vs_naive,
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
      values = c("skyblue", "gray", "tomato"), # to set the colours of our variable
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    coord_cartesian(ylim = c(0, max(-log10(df$pvalue))), xlim = c(-x_limit, x_limit)) +
    scale_x_continuous(breaks = seq(-x_limit, x_limit, 2)) +
    geom_text_repel(max.overlaps = Inf) # To show all labels
}


#' Generates a volcano plot using ggplot2
#' This function generates a base volcanoplot highlighting genes associated with a certain GOterm
#' that can then be expanded upon using further ggplot functions.
#'
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param res_enrich A enrichment result object created by for example using topGOtable
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param term_index The location (row) of your GO term of interest in your enrichment result
#' @param L2FC_cutoff A numeric value that sets the cutoff for the xintercept argument of ggplot
#' @param col_to_use The column in your differential expression results containing your gene symbols.
#'  If you don't have one it is created automatically
#' @param enrich_col column name from your res_enrich where the genes associated with your GOterm are stored (for example see the topGOtable result in mosdef)
#' @param down_col The colour for your downregulated genes, default is "gray"
#' @param up_col The colour for your upregulated genes, default is "gray"
#' @param highlight_col The colour for the genes associated with your GOterm default is "tomato"
#' @param overlaps number of overlaps ggrepel is supposed to allow when labeling
#'  (for more info check ggrepel documentation)
#'
#' @return A  \code{ggplot2} volcano plot object that can be extended upon by the user
#' @export
#'
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom AnnotationDbi mapIds
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline geom_point
#' theme_classic scale_color_manual coord_cartesian scale_x_continuous
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom rlang .data
#' @importFrom topGO runTest GenTable score sigGenes genesInTerm showSigOfNodes
#' annFUN annFUN.org
#'
#'
#'
#' @examples
#' library(ggplot2)
#' library(RColorBrewer) # for a colourful plot
#' library(ggrepel) # for nice annotations
#' library(DESeq2)
#' library("topGO")
#' library(AnnotationDbi)
#' library("org.Hs.eg.db")
#'
#'
#' data(res_de_macrophage, package = "mosdef")
#' data(res_enrich_macrophage_topGO, package = "mosdef")
#' p <- go_volcano(
#'   res_macrophage_IFNg_vs_naive,
#'   res_enrich = res_enrich_macrophage_topGO,
#'   term_index = 1,
#'   L2FC_cutoff = 1,
#'   mapping = "org.Hs.eg.db",
#'   overlaps = 20
#' )
#' p
go_volcano <- function(
    res_de,
    res_enrich,
    mapping = "org.Hs.eg.db",
    term_index,
    L2FC_cutoff = 1,
    col_to_use = NULL,
    enrich_col = "genes",
    down_col = "black",
    up_col = "black",
    highlight_col = "tomato",
    overlaps = 20) {
  if (is.null(col_to_use)) {
    res_de$symbol <- mapIds(get(mapping),
                            keys = row.names(res_de),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first"
    )
  } else {
    res_de$symbol <- res_de[[col_to_use]]
  }

  df <- deseqresult2df(res_de)




  # finding the highest value in the log2FoldChange column and rounding it up to get a nice symetric plot
  x_limit <- ceiling(max(abs(range(df$log2FoldChange, na.rm = TRUE))))



  df$diffexpressed <- "NO"

  df$diffexpressed[df$log2FoldChange > L2FC_cutoff & df$pvalue < 0.05] <- "UP"

  df$diffexpressed[df$log2FoldChange < -L2FC_cutoff & df$pvalue < 0.05] <- "DOWN"

  # calculate top 30 degenes based on pvalue
  # df$delabel <- ifelse(df$symbol %in% head(df[order(df$pvalue), "symbol"], labeled_genes), df$symbol, NA)


  # Either get indexes of the term name
  # grep("cell adhesion", res_enrich[["Term"]], fixed =TRUE)


  # or have the user tell us at which index the term of interest is

  test_vec <- res_enrich[[enrich_col]][term_index]
  test_vec <- strsplit(test_vec, ",")
  test_vec <- as.vector(test_vec)
  # test_for_plot <- test_vec[[1]][1:30]

  for (i in seq_len(length(df$id))) {
    if (df$symbol[i] %in% test_vec[[1]]) {
      df$de_label[i] <- df$symbol[i]
    } else {
      df$de_label[i] <- NA
    }
  }
  # df$delabel <- ifelse(df$symbol %in% head(df[order(df$de_label), "symbol"], labeled_genes), df$symbol, NA)

  p <- ggplot(data = df, aes(
    x = .data$log2FoldChange, y = -log10(.data$pvalue),
    colour = .data$diffexpressed, label = .data$de_label
  )) +
    geom_vline(xintercept = c(-L2FC_cutoff, L2FC_cutoff), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    geom_point() +
    theme_classic() +
    coord_cartesian(ylim = c(0, max(-log10(df$pvalue))), xlim = c(-x_limit, x_limit)) +
    scale_x_continuous(breaks = seq(-x_limit, x_limit, 2)) +
    geom_label_repel(max.overlaps = overlaps) # To show all labels
  # or   geom_text_repel(max.overlaps = Inf) not quite sure which is better


  # Define the genes to highlight
  genes_to_highlight <- c(test_vec[[1]])
  # Define a custom color scale
  custom_color_scale <- c("DOWN" = down_col, "NO" = "gray", "UP" = up_col, "Highlighted" = highlight_col)

  # Add another geom_point layer to highlight specific genes
  q <- p +
    geom_point(
      data = df[df$symbol %in% genes_to_highlight, ],
      aes(color = "Highlighted"), size = 3, shape = 16
    ) +
    scale_color_manual(
      values = custom_color_scale,
      labels = c("Downregulated", "GOterm", "Not significant", "Upregulated")
    )
}
