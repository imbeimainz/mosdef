#' Create a button for gene symbols in an Rmd
#'
#' A function to turn Gene Symbols into buttons in an Rmd linking to various Portals for further info
#' about these genes.
#' Current supported portals are: GeneCards, NCBI, GTEx, Uniprot, dbPTM, Human Protein Atlas
#'
#' @param df A dataframe with at least on column with gene Symbols named: SYMBOL
#' @param new_cols At least one of: "GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA" "PUBMED"
#' @param col_to_use name of the columns were the gene symbols are stored. Default is SYMBOL
#' @param output_format a parameter deciding which output format to return, either a DT:datatable (recommended)
#' or a simple dataframe (DF).In the latter case it is important that if the data is visualized with the
#'  \code{datatable} function the parameter escape must be set to FALSE
#' @param ens_col ame of the columns were the ensembl IDs are stored. 
#' @param ens_species The species you are working with to link to the correct gene on ensembl
#'  
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
#' data(res_de_macrophage, package = "mosdef")
#' res_de <- res_macrophage_IFNg_vs_naive
#' res_df <- deseqresult2df(res_de)
#' # Subsetting for quicker run
#' res_df <- res_df[1:100, ]
#' buttonifier(res_df)
#'
#' buttonifier(res_df,
#'   new_cols = c("NCBI", "HPA"),
#'   ens_col = "id",
#'   ens_species = "Homo_sapiens"
#' )
buttonifier <- function(df, new_cols = c("GC", "UNIPROT"),
                        col_to_use = "SYMBOL",
                        output_format = "DT",
                        ens_col = NULL,
                        ens_species = NULL) {
  if (!(col_to_use %in% colnames(df))) {
    stop(
      "The provided dataframe does not contain the column ", col_to_use, ". Please make ",
      "sure that there is a colum with gene symbols in your df and that its name is provided",
      " to the 'col_to_use' parameter. Please watch spelling as well as capital letters."
    )
  }



  if (!is.null(c(ens_col, ens_species))) {
    df[[ens_col]] <- create_link_ENS(df[[ens_col]], species = ens_species)
  } else if (!is.null(ens_col) & is.null(ens_species)) {
    warning(
      "Creating ensemble links requires an ID and the species you are analysing ",
      "You only provided an ID. "
    )
  } else if (is.null(ens_col) & !is.null(ens_species)) {
    warning(
      "Creating ensemble links requires an ID and the species you are analysing ",
      "You only provided a species."
    )
  }

  val <- df[[col_to_use]]

  match.arg(new_cols, choices = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED"), several.ok = TRUE)
  for (i in seq_len(length(new_cols))) {
    if ((new_cols[i] %in% c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED")) == FALSE) {
      warning(
        "Please make sure you used the values suggested in the documentation. \n",
        "One or more of the following values entered into new_cols is not supported: \n",
        new_cols[new_cols %in% c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA", "PUBMED") == FALSE]
      )
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

  if (output_format == "DT") {
    return(DT::datatable(df, escape = FALSE, rownames = FALSE))
  } else if (output_format == "DF") {
    df <- data.frame(df)
    return(df)
  }
}
