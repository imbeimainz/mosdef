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
#' @param overlaps number of overlaps ggrepel is supposed to allow when labelling
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
#'   res_enrich = topgoDE_macrophage_IFNg_vs_naive,
#'   term_index = 1,
#'   L2FC_cutoff = 1,
#'   mapping = "org.Hs.eg.db",
#'   overlaps = 30
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
    down_col ="black",
    up_col = "black",
    highlight_col = "tomato",
    overlaps = 20){
  
  if(is.null(col_to_use)){
    res_de$symbol <- mapIds(get(mapping),
                            keys = row.names(res_de),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first"
    )
  }else
    res_de$symbol <- res_de[[col_to_use]]
  
  df <- deseqresult2df(res_de)
  
  
  
  
  # finding the highest value in the log2FoldChange column and rounding it up to get a nice symetric plot
  x_limit <- ceiling(max(abs(range(df$log2FoldChange, na.rm = TRUE))))
  
  
  
  df$diffexpressed <- "NO"
  
  df$diffexpressed[df$log2FoldChange > L2FC_cutoff & df$pvalue < 0.05] <- "UP"
  
  df$diffexpressed[df$log2FoldChange < -L2FC_cutoff & df$pvalue < 0.05] <- "DOWN"
  
  # calculate top 30 degenes based on pvalue
  #df$delabel <- ifelse(df$symbol %in% head(df[order(df$pvalue), "symbol"], labeled_genes), df$symbol, NA)
  
  
  # Either get indexes of the term name
  #grep("cell adhesion", res_enrich[["Term"]], fixed =TRUE)
  
  
  #or have the user tell us at which index the term of interest is
  
  test_vec <- res_enrich[[enrich_col]][term_index]
  test_vec <- strsplit(test_vec, ",")
  test_vec <- as.vector(test_vec)
  #test_for_plot <- test_vec[[1]][1:30]
  
  for (i in 1:length(df$id)) {
    
    
    if (df$symbol[i] %in% test_vec[[1]]){
      df$de_label[i] <- df$symbol[i]
      
      
    } else{
      df$de_label[i] <- NA
    }
    
    
  }
  #df$delabel <- ifelse(df$symbol %in% head(df[order(df$de_label), "symbol"], labeled_genes), df$symbol, NA)
  
  p <-ggplot(data = df, aes(
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
    geom_point(data = df[df$symbol %in% genes_to_highlight, ],
               aes(color = "Highlighted"), size = 3, shape = 16) +
    scale_color_manual(values = custom_color_scale,
                       labels = c("Downregulated", "GOterm","Not significant", "Upregulated"))
  
}

