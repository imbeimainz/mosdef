#' Link to AMIGO database
#'
#' @param val Character, the GOID
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_GO("GO:0008150")
#'
create_link_GO <- function(val) {
  sprintf(
    '<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary">%s</a>',
    val,
    paste0(val, "@AMIGO")
  )
}

#' Link to the GeneCards database
#'
#' @param val Character, the gene symbol of interest
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_genecards("Oct4")
#'
create_link_genecards <- function(val) {
  sprintf(
    '<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    val, # link portion related to the gene
    .actionbutton_biocstyle, # button style
    paste0(val, "@GeneCards") # content of the button label
  )
}


#' Link to Pubmed
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_pubmed("Oct4")
#'
create_link_pubmed <- function(val) {
  paste0(
    '<a href="https://pubmed.ncbi.nlm.nih.gov/?term=', val, '" target="_blank" class="btn btn-primary">',
    paste0(val, "@Pubmed"), "</a>"
  )
}

#' Link to Ensemble database
#'
#' @param val Character, the gene symbol
#' @param species The species to be analyzed e.g "Mus_musculus"
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_ENS("ENSMUSG00000024406")
#'
create_link_ENS <- function(val, species = "Mus_musculus") {
  paste0('<a href="http://www.ensembl.org/', species, "/Gene/Summary?g=", val, '" target="_blank" class="btn btn-primary">', val, "</a>")
}

#' Link to NCBI database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_NCBI("Oct4")
#'
create_link_NCBI <- function(val) {
  paste0(
    '<a href="http://www.ncbi.nlm.nih.gov/gene/?term=', val, '[sym]" target="_blank" class="btn btn-primary">',
    paste0(val, "@NCBI"), "</a>"
  )
}

#' Link to the GTEx Portal
#'
#' @param val Character, the gene symbol of interest
#'
#' @return HTML for an action button
#' @export
#' @examples
#' create_link_GTEX("Oct4")
#'
create_link_GTEX <- function(val) {
  sprintf(
    '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
    val, # link portion related to the gene
    .actionbutton_biocstyle, # button style
    paste0(val, "@GTEX") # content of the button label
  )
}




#' Link to Uniprot database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_UniProt("Oct4")
#'
create_link_UniProt <- function(val) {
  sprintf(
    '<a href = "https://www.uniprot.org/uniprot/?query=%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
    val, # link portion related to the gene
    .actionbutton_biocstyle, # button style
    paste0(val, "@UNIPROT")
  ) # content of the button label
}

#' Link to dbPTM database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_dbPTM("Oct4")
#'
create_link_dbPTM <- function(val) {
  base_link_old <- "http://dbptm.mbc.nctu.edu.tw/"
  base_link_new <- "https://awi.cuhk.edu.cn/dbPTM/"


  sprintf(
    '<a href = "%s/info.php?id=%s_HUMAN" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-edit"></i>%s</a>',
    base_link_new, # main link to website
    val, # link portion related to the gene
    .actionbutton_biocstyle, # button style
    paste0(val, "@dbPTM")
  ) # content of the button label
}

#' Link to the Human Protein Atlas
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @export
#'
#' @examples
#' create_link_HPA("Oct4")
#'
create_link_HPA <- function(val) {
  sprintf(
    '<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
    val, # link portion related to the gene
    .actionbutton_biocstyle, # button style
    paste0(val, "@Human Protein Atlas")
  ) # content of the button label
}







#' Information on a GeneOntology identifier
#'
#' Assembles information, in HTML format, regarding a Gene Ontology identifier
#'
#' Also creates a link to the AmiGO database
#'
#' @param go_id Character, specifying the GeneOntology identifier for which
#' to retrieve information
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. If not provided, the experiment-related information is not
#' shown, and only some generic info on the identifier is displayed.
#'
#' @return HTML content related to a GeneOntology identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#'
#' @importFrom AnnotationDbi Term Ontology Definition Secondary  GOID
#' @importFrom htmltools tags HTML
#' @importFrom GO.db GOTERM
#'
#'
#' @examples
#' go_2_html("GO:0002250")
#' go_2_html("GO:0043368")
go_2_html <- function(go_id,
                      res_enrich = NULL) {
  fullinfo <- GO.db::GOTERM[[go_id]]
  if (is.null(fullinfo)) {
    return(HTML("Gene Ontology term not found!"))
  }
  # extracting the field/component values
  go_linkbutton <- create_link_GO(AnnotationDbi::GOID(fullinfo))
  go_term <- AnnotationDbi::Term(fullinfo)
  go_pubmed <- create_link_pubmed(go_term)
  go_ontology <- AnnotationDbi::Ontology(fullinfo)
  go_definition <- AnnotationDbi::Definition(fullinfo)
  go_synonims <- paste0(
    unlist(
      lapply(AnnotationDbi::Synonym(fullinfo), function(arg) {
        paste0(tags$b("Synonym: "), arg, tags$br())
      })
    ),
    collapse = ""
  )
  go_secondary <- Secondary(fullinfo)
  if (!is.null(res_enrich)) {
    go_pvalue <- res_enrich[(res_enrich$gs_id == go_id), "gs_pvalue"]
    go_zscore <- ifelse(
      "z_score" %in% colnames(res_enrich),
      format(round(res_enrich[(res_enrich$gs_id == go_id), "z_score"], 2), nsmall = 2),
      "NA - not yet computed"
    )
    go_aggrscore <- ifelse(
      "aggr_score" %in% colnames(res_enrich),
      format(round(res_enrich[(res_enrich$gs_id == go_id), "aggr_score"], 2), nsmall = 2),
      "NA - not yet computed"
    )
  }
  # assembling them together
  mycontent <- paste0(
    htmltools::tags$b("GO ID: "), go_linkbutton, htmltools::tags$br(),
    htmltools::tags$b("Pubmed results: "), go_pubmed, htmltools::tags$br(),
    htmltools::tags$b("Term: "), go_term, htmltools::tags$br(),
    ifelse(
      !is.null(res_enrich),
      paste0(htmltools::tags$b("p-value: "), go_pvalue, htmltools::tags$br(),
        htmltools::tags$b("Z-score: "), go_zscore, htmltools::tags$br(),
        htmltools::tags$b("Aggregated score: "), go_aggrscore, htmltools::tags$br(),
        collapse = ""
      ),
      ""
    ),
    tags$b("Ontology: "), go_ontology, tags$br(), htmltools::tags$br(),
    tags$b("Definition: "), go_definition, htmltools::tags$br(),
    go_synonims,
    ifelse(
      length(go_secondary) > 0,
      paste0(htmltools::tags$b("Secondary: "), go_secondary, collapse = ""),
      ""
    )
  )
  return(htmltools::HTML(mycontent))
}



#' Information on a gene
#'
#' Assembles information, in HTML format, regarding a gene symbol identifier
#'
#' Creates links to the NCBI and the GeneCards databases
#'
#' @param gene_id Character specifying the gene identifier for which to retrieve
#' information
#' @param res_de A `DESeqResults` object, storing the result of the differential
#' expression analysis. If not provided, the experiment-related information is not
#' shown, and only some generic info on the identifier is displayed.
#' The information about the gene is retrieved by matching on the `SYMBOL` column,
#' which should be provided in `res_de`.
#' @param col_to_use The column of your res_de object containing the gene symbols.
#' Default is "SYMBOL"
#'
#' @return HTML content related to a gene identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#'
#' @importFrom htmltools tags
#'
#'
#' @examples
#' geneinfo_2_html("ACTB")
#' geneinfo_2_html("Pf4")
geneinfo_2_html <- function(gene_id,
                            res_de = NULL,
                            col_to_use = "SYMBOL") {
  gene_ncbi_button <- create_link_NCBI(gene_id)
  gene_genecards_button <- create_link_genecards(gene_id)
  gene_gtex_button <- create_link_GTEX(gene_id)
  gene_uniProt_button <- create_link_UniProt(gene_id)
  gene_pubmed_button <- create_link_pubmed(gene_id)

  if (!is.null(res_de)) {
    gid <- match(gene_id, res_de[[col_to_use]])
    if (is.na(gid)) {
      message(
        "Could not find the specified gene (`", gene_id,
        "`) in the `res_de` object. \n",
        "Still, the general HTML content has been generated."
      )
      gene_adjpvalue <- htmltools::tags$em("not found")
      gene_logfc <- htmltools::tags$em("not found")
    } else {
      gene_adjpvalue <- format(res_de[gid, "padj"])
      gene_logfc <- format(round(res_de[gid, "log2FoldChange"], 2), nsmall = 2)
    }
  }

  mycontent <- paste0(
    htmltools::tags$b(gene_id), htmltools::tags$br(),
    "Link to the NCBI Gene database: ", gene_ncbi_button, htmltools::tags$br(),
    "Link to the GeneCards database: ", gene_genecards_button, htmltools::tags$br(),
    "Link to the GTEx Portal: ", gene_gtex_button, htmltools::tags$br(),
    "Link to the UniProt Portal: ", gene_uniProt_button, htmltools::tags$br(),
    "Link to related articles on Pubmed: ", gene_pubmed_button, htmltools::tags$br(),
    ifelse(
      !is.null(res_de),
      paste0(htmltools::tags$b("DE p-value (adjusted): "), gene_adjpvalue, htmltools::tags$br(),
        htmltools::tags$b("DE log2FoldChange: "), gene_logfc,
        collapse = ""
      ),
      ""
    )
  )
  return(htmltools::HTML(mycontent))
}



# Some constant values ----------------------------------------------------

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
