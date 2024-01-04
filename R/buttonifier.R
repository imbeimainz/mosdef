#' Create a button for genesymbols in an RMD
#'
#' A function to turn Gene Symbols into buttons in an RMD linking to various Portals for further info
#' about these genes.
#' Current supported portals are: Genecards, NCBI, GTEX, Uniprot, dbPTM, Human Protein Atlas
#'
#' @param df A dataframe with at least on coloumn with gene Symbols named: SYMBOL
#' @param new_cols At least one of: "GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"
#' @param col_to_use name of the coloumns were the gene symbols are stored. Default is SYMBOL
#' @param output_format a parameter deciding which output format to return, either a DT:datatable (recommended)
#' or a simple dataframe (DF).In the latter case it is important that if the data is visualized with the
#'  \code{datatable} function the parameter escape must be set to FALSE
#'
#' @return A dataframe or a \code{DT} datatable object with columns adding HTML objects that link to websites with further information on the genes in question.
#' @export
#' 
#' 
#' 
#' @importFrom DT datatable
#'
#' @examples
#' library(dplyr)
#' library(DESeq2)
#' data("gse", package = "macrophage")
#'
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
# changing the ids to Ensembl instead of the Gencode used in the object
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#' dds_macrophage
#' dds_macrophage <- DESeq(dds_macrophage)
#'
#' res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
#'                                         contrast = c("condition", "IFNg", "naive"),
#'                                         lfcThreshold = 1, alpha = 0.05)
#' res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL
#' res_df <- as.data.frame(res_macrophage_IFNg_vs_naive@listData)
#' res_df <-res_df[1:100,]
#' buttonifier(res_df)
buttonifier <- function(df, new_cols = c("GC", "UNIPROT"), col_to_use = "SYMBOL", output_format = "DT"){
  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
  val <- df[[col_to_use]]
  match.arg(new_cols,choices = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"), several.ok = TRUE)
  match.arg(output_format, choices = c("DT", "DF"))

  if( "GC" %in% new_cols){
    df$SYMBOL_GC <- create_link_genecards(df[[col_to_use]])

  }
  if( "NCBI" %in% new_cols){
    df$SYMBOL_NCBI <- create_link_NCBI(df[[col_to_use]])

  }

  if( "GTEX" %in% new_cols){
    df$SYMBOL_GTEX <- create_link_GTEX(df[[col_to_use]])

  }

  if( "UNIPROT" %in% new_cols){
    df$SYMBOL_UNIPROT <- create_link_UniProt(df[[col_to_use]])

  }

  if( "dbPTM" %in% new_cols){
    df$SYMBOL_dbPTM <- create_link_dbPTM(df[[col_to_use]])

  }

  if( "HPA" %in% new_cols){
    df$SYMBOL_HPA <- create_link_HPA(df[[col_to_use]])

  }

  df <- df %>%
    select(-SYMBOL)
  if (output_format == "DT"){
    return(DT::datatable(df, escape = FALSE))

  }else if (output_format == "DF"){

    return(df)
  }





}


















