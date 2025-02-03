#' Get an annotation data frame from org db packages
#'
#' @param de_container An object containing the data for a Differential
#' Expression workflow (e.g. `DESeq2`, `edgeR` or `limma`).
#' Currently, this can be a `DESeqDataSet` object, normally obtained after
#' running your data through the `DESeq2` framework.
#' @param orgdb_package Character string, named as the `org.XX.eg.db`
#' package which should be available in Bioconductor
#' @param id_type Character, the ID type of the genes as in the row names of
#' the `de_container`, to be used in the call to [AnnotationDbi::mapIds()]
#' @param key_for_genenames Character, corresponding to the column name for the
#' key in the orgDb package containing the official gene name (often called
#' gene symbol).
#' This parameter defaults to "SYMBOL", but can be adjusted in case the key is not
#' found in the annotation package (e.g. for `org.Sc.sgd.db`).
#'
#' @return A data frame to be used for annotation of genes, with the main
#' information encoded in the `gene_id` and `gene_name` columns.
#'
#' @importFrom AnnotationDbi mapIds keytypes
#'
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#'
#' # dds object
#' data(gse, package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#'
#' anno_df <- get_annotation_orgdb(dds_macrophage, "org.Hs.eg.db", "ENSEMBL")
#'
#' head(anno_df)
get_annotation_orgdb <- function(de_container,
                                 orgdb_package,
                                 id_type,
                                 key_for_genenames = "SYMBOL") {

  ## TODO: see if we need any checks
  if (is(de_container, "DESeqDataSet")) {
    ## TODO
    # all we need is rownames
  }



  if (is.null(orgdb_package))
    stop("Select an annotation package to generate the corresponding annotation")

  orgdbpkgs <- data.frame(
    pkg = c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
            "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Hs.ipi.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db",
            "org.Pt.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db", "org.Xl.eg.db"),
    descr = c("Anopheles", "Arabidopsis", "Bovine", "Worm", "Canine", "Fly", "Zebrafish", "E coli strain K12", "E coli strain Sakai", "Chicken",
              "Human", "org.Hs.ipi.db", "Mouse", "Rhesus", "Malaria", "Chimp", "Rat", "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
              "Xenopus"),
    stringsAsFactors = FALSE
  )

  if (!(orgdb_package %in% orgdbpkgs$pkg)) {
    message("The orgDB package is most likely not existent in Bioconductor")
    message("It should be one of", orgdbpkgs$pkg)
  }

  if (!require(orgdb_package, character.only = TRUE))
    stop("The package ", orgdb_package, " is not installed/available. Try installing it with BiocManager::install('", orgdb_package, "')")

  if (!(id_type %in% keytypes(eval(parse(text = orgdb_package))))) {
    stop("The key you provided is not listed as key for the annotation package. Please try one of ",
         paste(keytypes(eval(parse(text = orgdb_package))), collapse = ","))
  }

  if (!(key_for_genenames %in% keytypes(eval(parse(text = orgdb_package))))) {
    stop("The key specified for containing gene names is not included in the annotation package. Please try one of ",
         paste(keytypes(eval(parse(text = orgdb_package))), collapse = ","))
  }

  pkg <- eval(parse(text = orgdb_package))

  if (id_type == "SYMBOL")
    warning("You probably do not need to convert symbol to symbol") # the performance would somehow be affected

  anns_vec <- mapIds(pkg,
                     keys = rownames(de_container),
                     column = key_for_genenames,
                     keytype = id_type)

  anno_df <- data.frame(
    gene_id = rownames(de_container),
    gene_name = anns_vec,
    stringsAsFactors = FALSE,
    row.names = rownames(de_container)
  )
  return(anno_df)
}
