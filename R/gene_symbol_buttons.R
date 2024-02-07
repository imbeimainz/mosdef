#' Create a button for gene symbols in an RMD
#'
#' A function to turn Gene Symbols into buttons in an RMD linking to various Portals for further info
#' about these genes.
#' Current supported portals are: Genecards, NCBI, GTEX, Uniprot, dbPTM, Human Protein Atlas
#'
#' @param df A dataframe with at least on column with gene Symbols named: SYMBOL
#' @param new_cols At least one of: "GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"
#' @param col_to_use name of the columns were the gene symbols are stored. Default is SYMBOL
#' @param output_format a parameter deciding which output format to return, either a DT:datatable (recommended)
#' or a simple dataframe (DF).In the latter case it is important that if the data is visualized with the
#'  \code{datatable} function the parameter escape must be set to FALSE
#'
#' @return TODO
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' library(dplyr)
#' library(DESeq2)
#' data(res_airway, package = "mosdef")
#' res_de <- res_airway
#' res_df <- deseqresult2df(res_de)
#' # Subsetting for quicker run
#' res_df <- res_df[1:100,]
#' gene_symbol_buttons(res_df)
gene_symbol_buttons <- function(df, new_cols = c("GC", "UNIPROT"), col_to_use = "SYMBOL", output_format = "DT") {
  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
  val <- df[[col_to_use]]
  match.arg(new_cols, choices = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"), several.ok = TRUE)
  match.arg(output_format, choices = c("DT", "DF"))
  # GeneCards
  if ("GC" %in% new_cols) {
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_GC[i] <- sprintf(
        '<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn     btn-primary" style = "%s">%s</a>',
        val[i],
        .actionbutton_biocstyle,
        paste0(val[i], "@GeneCards")
      )
    }
  }

  # NCBI
  if ("NCBI" %in% new_cols) {
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_NCBI[i] <- paste0(
        '<a href="http://www.ncbi.nlm.nih.gov/gene/?term=', val[i], '[sym]" target="_blank" class="btn btn-primary">',
        paste0(val[i], "@NCBI"), "</a>"
      )
    }
  }

  # GTEx
  if ("GTEX" %in% new_cols) {
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_GTEX[i] <- sprintf(
        '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
        val[i],
        .actionbutton_biocstyle,
        paste0(val[i], "@GTEX")
      )
    }
  }

  # Uniprot
  if ("UNIPROT" %in% new_cols) {
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_UniProt[i] <- sprintf(
        '<a href = "https://www.uniprot.org/uniprot/?query=%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
        val[i],
        .actionbutton_biocstyle,
        paste0(val[i], "@UNIPROT")
      )
    }
  }

  # dbPTM
  if ("dbPTM" %in% new_cols) {
    base_link_old <- "http://dbptm.mbc.nctu.edu.tw/"
    base_link_new <- "https://awi.cuhk.edu.cn/dbPTM/"

    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_dbPTM[i] <- sprintf(
        '<a href = "%s/info.php?id=%s_HUMAN" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-edit"></i>%s</a>',
        base_link_new, # main link to website
        val[i], # link portion related to the gene
        .actionbutton_biocstyle, # button style
        paste0(val[i], "@dbPTM")
      ) # content of the button label
    }
  }

  # Human protein atlas Portal
  if ("HPA" %in% new_cols) {
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_HPA[i] <- sprintf(
        '<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
        val[i],
        .actionbutton_biocstyle,
        paste0(val[i], "@Human Protein Atlas")
      )
    }
  }
  df <- df |>
    select(-.data$SYMBOL)
  if (output_format == "DT") {
    return(DT::datatable(df, escape = FALSE))
  } else if (output_format == "DF") {
    return(df)
  }
}
