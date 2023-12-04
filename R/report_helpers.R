create_link_GO <- function(val) {
  sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary">%s</a>',val,val)
}

create_link_ENS  <- function(val, species="Mus_musculus") {
  paste0('<a href="http://www.ensembl.org/',species,'/Gene/Summary?g=',val,'" target="_blank" class="btn btn-primary">',val,'</a>')
}

#' Link to NCBI database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
create_link_NCBI <- function(val) {
  # possibilities:
  # ncbi
  # genecards
  paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',val,'[sym]" target="_blank" class="btn btn-primary">',val,'</a>')
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
#' See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#'
#' @return HTML content related to a GeneOntology identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' go_2_html("GO:0002250")
#' go_2_html("GO:0043368")
go_2_html <- function(go_id,
                      res_enrich = NULL) {
  fullinfo <- GOTERM[[go_id]]
  if (is.null(fullinfo)) {
    return(HTML("Gene Ontology term not found!"))
  }
  # extracting the field/component values
  go_linkbutton <- .link2amigo(GOID(fullinfo))
  go_term <- Term(fullinfo)
  go_ontology <- Ontology(fullinfo)
  go_definition <- Definition(fullinfo)
  go_synonims <- paste0(
    unlist(
      lapply(Synonym(fullinfo), function(arg) {
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
    tags$b("GO ID: "), go_linkbutton, tags$br(),
    tags$b("Term: "), go_term, tags$br(),
    ifelse(
      !is.null(res_enrich),
      paste0(tags$b("p-value: "), go_pvalue, tags$br(),
             tags$b("Z-score: "), go_zscore, tags$br(),
             tags$b("Aggregated score: "), go_aggrscore, tags$br(),
             collapse = ""
      ),
      ""
    ),
    tags$b("Ontology: "), go_ontology, tags$br(), tags$br(),
    tags$b("Definition: "), go_definition, tags$br(),
    go_synonims,
    ifelse(
      length(go_secondary) > 0,
      paste0(tags$b("Secondary: "), go_secondary, collapse = ""),
      ""
    )
  )
  return(HTML(mycontent))
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
#'
#' @return HTML content related to a gene identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' geneinfo_2_html("ACTB")
#' geneinfo_2_html("Pf4")
geneinfo_2_html <- function(gene_id,
                            res_de = NULL) {
  gene_ncbi_button <- .link2ncbi(gene_id)
  gene_genecards_button <- .link2genecards(gene_id)
  gene_gtex_button <- .link2gtex(gene_id)

  if (!is.null(res_de)) {
    gid <- match(gene_id, res_de$SYMBOL)
    if (is.na(gid)) {
      message(
        "Could not find the specified gene (`", gene_id,
        "`) in the `res_de` object. \n",
        "Still, the general HTML content has been generated."
      )
      gene_adjpvalue <- tags$em("not found")
      gene_logfc <- tags$em("not found")
    } else {
      gene_adjpvalue <- format(res_de[gid, "padj"])
      gene_logfc <- format(round(res_de[gid, "log2FoldChange"], 2), nsmall = 2)
    }
  }

  mycontent <- paste0(
    tags$b(gene_id), tags$br(),
    "Link to the NCBI Gene database: ", gene_ncbi_button, tags$br(),
    "Link to the GeneCards database: ", gene_genecards_button, tags$br(),
    "Link to the GTEx Portal: ", gene_gtex_button, tags$br(),
    ifelse(
      !is.null(res_de),
      paste0(tags$b("DE p-value (adjusted): "), gene_adjpvalue, tags$br(),
             tags$b("DE log2FoldChange: "), gene_logfc,
             collapse = ""
      ),
      ""
    )
  )
  return(HTML(mycontent))
}



#' Link to the GeneCards database
#'
#' @param val Character, the gene symbol of interest
#'
#' @return HTML for an action button
#' @noRd
create_link_GC <- function(val) {
  sprintf(
    '<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}

#' Link to the GTEx Portal
#'
#' @param val Character, the gene symbol of interest
#'
#' @return HTML for an action button
#' @noRd
create_link_GTEX <- function(val) {
  sprintf(
    '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}







#' Link to the Uniprot Portal
create_link_UniProt <- function(val) {
  sprintf('<a href = "https://www.uniprot.org/uniprot/?query=%s&sort=score" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the dbPTM Portal
create_link_dbPTM <- function(val) {
  sprintf('<a href = "http://dbptm.mbc.nctu.edu.tw/info.php?id=%s_HUMAN" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-edit"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the human protein atlas Portal
create_link_HPA <- function(val) {
  sprintf('<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
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
#'
#' @return HTML content related to a gene identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' geneinfo_2_html("ACTB")
#' geneinfo_2_html("Pf4")
geneinfo_2_html <- function(gene_id,
                            res_de = NULL) {
  gene_ncbi_button <- .link2ncbi(gene_id)
  gene_genecards_button <- .link2genecards(gene_id)
  gene_gtex_button <- .link2gtex(gene_id)
  
  if (!is.null(res_de)) {
    gid <- match(gene_id, res_de$SYMBOL)
    if (is.na(gid)) {
      message(
        "Could not find the specified gene (`", gene_id,
        "`) in the `res_de` object. \n",
        "Still, the general HTML content has been generated."
      )
      gene_adjpvalue <- tags$em("not found")
      gene_logfc <- tags$em("not found")
    } else {
      gene_adjpvalue <- format(res_de[gid, "padj"])
      gene_logfc <- format(round(res_de[gid, "log2FoldChange"], 2), nsmall = 2)
    }
  }
  
  mycontent <- paste0(
    tags$b(gene_id), tags$br(),
    "Link to the NCBI Gene database: ", gene_ncbi_button, tags$br(),
    "Link to the GeneCards database: ", gene_genecards_button, tags$br(),
    "Link to the GTEx Portal: ", gene_gtex_button, tags$br(),
    ifelse(
      !is.null(res_de),
      paste0(tags$b("DE p-value (adjusted): "), gene_adjpvalue, tags$br(),
             tags$b("DE log2FoldChange: "), gene_logfc,
             collapse = ""
      ),
      ""
    )
  )
  return(HTML(mycontent))
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
#' See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#'
#' @return HTML content related to a GeneOntology identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' go_2_html("GO:0002250")
#' go_2_html("GO:0043368")
go_2_html <- function(go_id,
                      res_enrich = NULL) {
  fullinfo <- GOTERM[[go_id]]
  if (is.null(fullinfo)) {
    return(HTML("Gene Ontology term not found!"))
  }
  # extracting the field/component values
  go_linkbutton <- .link2amigo(GOID(fullinfo))
  go_term <- Term(fullinfo)
  go_ontology <- Ontology(fullinfo)
  go_definition <- Definition(fullinfo)
  go_synonims <- paste0(
    unlist(
      lapply(Synonym(fullinfo), function(arg) {
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
    tags$b("GO ID: "), go_linkbutton, tags$br(),
    tags$b("Term: "), go_term, tags$br(),
    ifelse(
      !is.null(res_enrich),
      paste0(tags$b("p-value: "), go_pvalue, tags$br(),
             tags$b("Z-score: "), go_zscore, tags$br(),
             tags$b("Aggregated score: "), go_aggrscore, tags$br(),
             collapse = ""
      ),
      ""
    ),
    tags$b("Ontology: "), go_ontology, tags$br(), tags$br(),
    tags$b("Definition: "), go_definition, tags$br(),
    go_synonims,
    ifelse(
      length(go_secondary) > 0,
      paste0(tags$b("Secondary: "), go_secondary, collapse = ""),
      ""
    )
  )
  return(HTML(mycontent))
}
# Some constant values ----------------------------------------------------

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"


