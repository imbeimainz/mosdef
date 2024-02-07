#' Create a button for gene symbols in an RMD
#'
#' A function to turn Gene Symbols into buttons in an RMD linking to various Portals for further info
#' about these genes.
#' Current supported portals are: Genecards, NCBI, GTEX, Uniprot, dbPTM, Human Protein Atlas
#'
#' @param df A dataframe with at least on column with gene Symbols named: SYMBOL
#' @param new_cols At least one of: "GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA" "PUBMED"
#' @param col_to_use name of the columns were the gene symbols are stored. Default is SYMBOL
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
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
#' @examples
#' data(res_airway, package = "mosdef")
#' res_de <- res_airway
#' res_df <- deseqresult2df(res_de)
#' # Subsetting for quicker run
#' res_df <- res_df[1:100,]
#' buttonifier(res_df)
buttonifier <- function(df, new_cols = c("GC", "UNIPROT"), col_to_use = "SYMBOL", output_format = "DT") {
  if (!(col_to_use %in% colnames(df))) {
    stop(
      "The provided dataframe does not contain the column ", col_to_use, ". Please make ",
      "sure that ther is a colum with gene symbols in your df and that its name is provided",
      " to the 'col_to_use' parameter. Please watch spelling as well as capital letters."
    )
  }


  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
  val <- df[[col_to_use]]

  match.arg(new_cols, choices = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED"), several.ok = TRUE)
  for (i in 1:length(new_cols)) {
    if ((new_cols[i] %in% c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED")) == FALSE) {
      warning(paste0(
        "Please make sure you used the values suggested in the documentation. \n",
        "One or more of the following values entered into new_cols is not supported: \n",
        new_cols[new_cols %in% c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED") == FALSE]
      ))
    }
  }


  match.arg(output_format, choices = c("DT", "DF"))

  if ("GC" %in% new_cols) {
    df$SYMBOL_GC <- create_link_genecards(df[[col_to_use]])
  }
  if ("NCBI" %in% new_cols) {
    df$SYMBOL_NCBI <- create_link_NCBI(df[[col_to_use]])
  }

  if ("GTEX" %in% new_cols) {
    df$SYMBOL_GTEX <- create_link_GTEX(df[[col_to_use]])
  }

  if ("UNIPROT" %in% new_cols) {
    df$SYMBOL_UNIPROT <- create_link_UniProt(df[[col_to_use]])
  }

  if ("dbPTM" %in% new_cols) {
    df$SYMBOL_dbPTM <- create_link_dbPTM(df[[col_to_use]])
  }

  if ("HPA" %in% new_cols) {
    df$SYMBOL_HPA <- create_link_HPA(df[[col_to_use]])
  }
  
  if ("PUBMED" %in% new_cols) {
    df$SYMBOL_PUBM <- create_link_pubmed(df[[col_to_use]])
  }

  df <- df |>
    select(-.data[[col_to_use]])
  if (output_format == "DT") {
    return(DT::datatable(df, escape = FALSE))
  } else if (output_format == "DF") {
    df <- data.frame(df)
    return(df)
  }
}
