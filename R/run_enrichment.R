#' Extract functional terms enriched in the DE genes, based on topGO
#'
#' A wrapper for extracting functional GO terms enriched in the DE genes, based on
#' the algorithm and the implementation in the topGO package
#'
#' Allowed values assumed by the \code{topGO_method2} parameter are one of the
#' following: \code{elim}, \code{weight}, \code{weight01}, \code{lea},
#' \code{parentchild}. For more details on this, please refer to the original
#' documentation of the \code{topGO} package itself
#'
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param dds A DESeqDataset object created using \code{DESeq2}
#' @param de_genes A vector of (differentially expressed) genes
#' @param bg_genes A vector of background genes, e.g. all (expressed) genes in the assays
#' @param top_de numeric, how many of the top differentially expressed genes to use for the enrichment analysis.
#'  Attempts to reduce redundancy. Assumes the data is sorted by padj (default in DESeq2).
#' @param min_counts numeric, min number of counts a gene needs to have to be included
#' in the geneset that the de genes are compared to. Default is 0, recommended only for advanced users.
#' @param ontology Which Gene Ontology domain to analyze: \code{BP} (Biological Process), \code{MF} (Molecular Function), or \code{CC} (Cellular Component)
#' @param annot Which function to use for annotating genes to GO terms. Defaults to \code{annFUN.org}
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param geneID Which format the genes are provided. Defaults to \code{symbol}, could also be
#' \code{entrez} or \code{ENSEMBL}
#' @param full_names_in_rows Logical, whether to display or not the full names for the GO terms
#' @param de_type One of: 'up', 'down', or 'up_and_down' Which genes to use for GOterm calculations:
#'  upregulated, downregulated or both
#' @param add_gene_to_terms Logical, whether to add a column with all genes annotated to each GO term
#' @param output_file Name of the file the result should be written into
#' @param topGO_method2 Character, specifying which of the methods implemented by \code{topGO} should be used, in addition to the \code{classic} algorithm. Defaults to \code{elim}
#' @param do_padj Logical, whether to perform the adjustment on the p-values from the specific
#' topGO method, based on the FDR correction. Defaults to FALSE, since the assumption of
#' independent hypotheses is somewhat violated by the intrinsic DAG-structure of the Gene
#' Ontology Terms
#' @param verbose Logical, whether to add messages telling the user which steps were taken
#'
#' @importFrom methods is new
#' @importFrom AnnotationDbi mapIds Term
#' @importFrom topGO runTest GenTable score sigGenes genesInTerm showSigOfNodes
#' annFUN annFUN.org
#' @importFrom utils write.table
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#'
#' @family Enrichment functions
#'
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' data(gse, package = "macrophage")
#'
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#' dds_macrophage <- DESeq(dds_macrophage)
#'
#' data(res_de_macrophage, package = "mosdef")
#'
#' library("AnnotationDbi")
#' library("org.Hs.eg.db")
#' library("topGO")
#' topgoDE_macrophage <- topGOtable(
#'   res_de = res_macrophage_IFNg_vs_naive,
#'   dds = dds_macrophage,
#'   ontology = "BP",
#'   mapping = "org.Hs.eg.db",
#'   geneID = "symbol",
#' )
run_topGO <- function(res_de = NULL, # Differentially expressed genes
                      dds = NULL, # background genes
                      de_genes = NULL,
                      bg_genes = NULL,
                      top_de = NULL,
                      min_counts = 0,
                      ontology = "BP", # could use also "MF"
                      annot = annFUN.org, # parameters for creating topGO object
                      mapping = "org.Mm.eg.db",
                      geneID = "symbol", # could also beID = "entrez"
                      # topTablerows = 200,     was used to limit table size. Now just use nrows
                      full_names_in_rows = TRUE,
                      add_gene_to_terms = TRUE,
                      de_type = "up_and_down",
                      # plot_graph = FALSE,
                      # plot_nodes = 10,
                      # write_output = FALSE,
                      output_file = "",
                      topGO_method2 = "elim",
                      do_padj = FALSE,
                      verbose = TRUE) {
  ## Checks:

  # Check if de-type is correct
  # if(!(de_type %in% c("up_and_down","up", "down")))

  match.arg(de_type, choices = c("up_and_down", "up", "down"), several.ok = FALSE)

  # stop("The de_type argument must be one of: 'up_and_down', 'up', 'down'")


  # Check if there is any input at all
  if (is.null(c(de_genes, bg_genes, dds, res_de))) {
    stop(
      "Please provide one of the following forms of input: \n",
      "A vector of differentially expressed genes and a vector of backgoud genes. \n",
      "A DEseqDataset  and a DESeqResults generated by DESeq2."
    )
  }

  # Check if there only a res_de is given
  if (!is.null(res_de) & is.null(dds)) {
    stop("Please also provide a DESeq2Dataset (dds) object.")
  }

  # check if only dds is given
  if (!is.null(dds) & is.null(res_de)) {
    stop("Please also provide a DESeq2 result object.")
  }

  # Check if there only a bg but no de genes is given
  if (!is.null(bg_genes) & is.null(de_genes)) {
    stop("Please also provide a vector of  differentially expressed genes.")
  }

  # Check if there is only a de but no bg genes given
  if (!is.null(de_genes) & is.null(bg_genes)) {
    stop("Please also provide  a vector of background genes.")
  }


  # checking the additional topGO_method2
  topgo_methods <- c("elim", "weight", "weight01", "lea", "parentchild")
  if (!(topGO_method2 %in% topgo_methods)) {
    stop(
      "Please provide one of the following topGO methods in addition to the classic method:\n",
      paste(topgo_methods, collapse = ", ")
    )
  }

  if ((de_type == "up" | de_type == "down") && !is.null(de_genes)) {
    stop(
      "The argument de_type can only be used if a dds and a res_de object are provided:\n",
      "please either provide these objects or if you want to work with gene vectors set de_type to: 'up_and_down'"
    )
  }

  annot_to_map_to <- get(mapping)

  if (!is.null(res_de) && !is.null(dds)) {
    # Check if the inputs are the correct type

    if (!is(dds, "DESeqDataSet")) {
      stop("The provided `dds` is not a DESeqDataSet object, please check your input parameters.")
    }

    if (!is(res_de, "DESeqResults")) {
      stop("The provided `res_de` is not a DESeqResults object, please check your input parameters.")
    }

    # checking that results and dds are related
    ## at least a subset of dds should be in res
    if (!all(rownames(res_de) %in% rownames(dds))) {
      warning(
        "It is likely that the provided `dds` and `res_de` objects are not related ",
        "to the same dataset (the row names of the results are not all in the dds). ",
        "Are you sure you want to proceed?"
      )
    }

    # Check if DESeq was run on the dds
    if (!"results" %in% mcols(mcols(dds))$type) {
      stop("I couldn't find results in your dds. You should first run DESeq2::DESeq() on your dds.")
    }

    res_de$symbol <- AnnotationDbi::mapIds(annot_to_map_to,
                                           keys = row.names(res_de),
                                           column = "SYMBOL",
                                           keytype = "ENSEMBL",
                                           multiVals = "first"
    )

    resOrdered <- as.data.frame(res_de[order(res_de$padj), ])

    if (de_type == "up_and_down") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
    } else if (de_type == "up") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
      de_df <- de_df[de_df$log2FoldChange >= 0, ]
    } else if (de_type == "down") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
      de_df <- de_df[de_df$log2FoldChange <= 0, ]
    }
    de_genes <- de_df$symbol
    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]
    }
    bg_ids <- rownames(dds)[rowSums(counts(dds)) > min_counts]
    bg_genes <- AnnotationDbi::mapIds(annot_to_map_to,
                                      keys = bg_ids,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first"
    )
    if (verbose) {
      message(
        "Your dataset has ",
        nrow(de_df),
        " DE genes. You selected ",
        length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / nrow(de_df)) * 100), # sprintf format with 2 decimal places
        ") genes. You analysed all ",
        de_type,
        "-regulated genes"
      )
    }
  } else if (!is.null(c(bg_genes, de_genes))) {
    all_de <- length(de_genes)

    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]
    }
    if (verbose) {
      message(
        "Your dataset has ",
        all_de,
        " DE genes.You selected ",
        length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / all_de) * 100), # sprintf format with 2 decimal places
        ") genes. You analysed all ",
        de_type,
        "-regulated genes"
      )
    }
  }

  # creating the vectors
  de_genes_input <- factor(as.integer(bg_genes %in% de_genes))
  names(de_genes_input) <- bg_genes

  # instantiating the topGOdata object
  suppressMessages({
    GOdata <- new("topGOdata",
                  ontology = ontology,
                  allGenes = de_genes_input,
                  nodeSize = 10,
                  annot = annot,
                  mapping = mapping,
                  ID = geneID
    )
  })
  # performing the test(s)
  suppressMessages({
    result_method2 <- runTest(GOdata, algorithm = topGO_method2, statistic = "fisher")
  })
  suppressMessages({
    resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  })
  sTab <- GenTable(GOdata,
                   p.value_method2 = result_method2,
                   p.value_classic = resultClassic,
                   orderBy = "p.value_method2",
                   ranksOf = "p.value_classic",
                   topNodes = length(score(resultClassic))
  )

  names(sTab)[which(names(sTab) == "p.value_method2")] <- paste0("p.value_", topGO_method2)

  sTab[["p.value_classic"]] <- as.numeric(sTab[["p.value_classic"]])
  sTab[[paste0("p.value_", topGO_method2)]] <- as.numeric(sTab[[paste0("p.value_", topGO_method2)]])

  # if FDR, then apply it here
  if (do_padj) {
    sTab[[paste0("padj_BY_", topGO_method2)]] <- p.adjust(sTab[[paste0("p.value_", topGO_method2)]], method = "BY")
  }

  # subset to specified number of rows
  ### topTablerows <- min(nrow(sTab), topTablerows)
  topTablerows <- nrow(sTab)
  sTab <- sTab[seq_len(topTablerows), ]

  if (full_names_in_rows) {
    sTab$Term <- vapply(sTab$GO.ID, function(go) {
      Term(GOTERM[[go]])
    }, FUN.VALUE = character(1))
  }

  if (add_gene_to_terms) {
    # adapted from an elegant one liner found here: https://support.bioconductor.org/p/65856/
    SignificantGenes <- sigGenes(GOdata)
    sTab$genes <- lapply(sTab$GO.ID, function(x) {
      genes <- genesInTerm(GOdata, x)
      tmp <- genes[[1]][genes[[1]] %in% SignificantGenes]
    })
    # coerce the list to a comma separated vector
    sTab$genes <- unlist(lapply(sTab$genes, function(arg) paste(arg, collapse = ",")))
  }
  # message for filters or a summary here
  message(
    nrow(sTab),
    " GO terms were analyzed. Not all of them are significantly enriched.\n",
    "We suggest further subsetting the output list by for example: \n",
    "using a pvalue cutoff in the column: \n",
    "'p.value_elim'."
  )

  return(sTab)
}


#' Extract functional terms enriched in the DE genes, based on goseq
#'
#' A wrapper for extracting functional GO terms enriched in a list of (DE) genes,
#' based on the algorithm and the implementation in the goseq package
#'
#' Note: the feature length retrieval is based on the \code{\link{goseq}} function,
#' and requires that the corresponding TxDb packages are installed and available
#'
#' @param de_genes A vector of (differentially expressed) genes
#' @param bg_genes A vector of background genes, e.g. all (expressed) genes
#' in the assays
#' @param dds A DESeqDataset object created using \code{DESeq2}
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param top_de numeric, how many of the top differentially expressed genes to use for the enrichment analysis.
#'  Attempts to reduce redundancy. Assumes the data is sorted by padj (default in DESeq2).
#' @param min_counts numeric, min number of counts a gene needs to have to be included
#' in the geneset that the de genes are compared to. Default is 0, recommended only for advanced users.
#' @param genome A string identifying the genome that genes refer to, as in the
#' \code{\link{goseq}} function
#' @param id A string identifying the gene identifier used by genes, as in the
#' \code{\link{goseq}} function
#' @param de_type One of: 'up', 'down', or 'up_and_down' Which genes to use for GOterm calculations:
#'  upregulated, downregulated or both
#' @param testCats A vector specifying which categories to test for over representation amongst DE genes - can be any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG"
#' @param FDR_GO_cutoff Numeric value for subsetting the results
#' @param mapping Character string, named as the \code{org.XX.eg.db}
#' package which should be available in Bioconductor
#' @param add_gene_to_terms Logical, whether to add a column with all genes annotated
#' to each GO term
#' @param verbose Logical, whether to add messages telling the user which steps were taken
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#'
#' @export
#'
#' @importFrom stats p.adjust
#' @importFrom goseq nullp goseq getgo
#' @importFrom AnnotationDbi mapIds
#'
#' @family Enrichment functions
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' data(gse, package = "macrophage")
#'
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#' dds_macrophage <- DESeq(dds_macrophage)
#'
#' data(res_de_macrophage, package = "mosdef")
#' res_de <- res_macrophage_IFNg_vs_naive
#' mygo <- goseqTable(
#'   res_de = res_macrophage_IFNg_vs_naive,
#'   dds = dds_macrophage,
#'   mapping = "org.Hs.eg.db",
#'   testCats = "GO:BP",
#'   add_gene_to_terms = TRUE
#' )
#'
#' head(mygo)
run_goseq <- function(res_de = NULL,
                      dds = NULL,
                      de_genes = NULL, # Differentially expressed genes
                      bg_genes = NULL, # background genes, normally = rownames(cds) or filtering to genes,
                      top_de = NULL,
                      min_counts = 0,
                      #  with at least 1 read - could also be ls(org.Mm.egGO)
                      genome = "hg38",
                      id = "ensGene",
                      de_type = "up_and_down",
                      testCats = c("GO:BP", "GO:MF", "GO:CC"),
                      FDR_GO_cutoff = 1,
                      # nTop = 200,
                      mapping = "org.Hs.eg.db",
                      # testKegg=TRUE,
                      # keggObject=mapPathwayToName("mmu"), # need the dedicated function!!
                      # writeOutput=FALSE,
                      add_gene_to_terms = TRUE,
                      verbose = TRUE
                      # outputFiles_goseq="",outputFiles_goseq_kegg=""
                      ## TODO TODO: bring back in action the function
                      ## add genes annotated to each term
                      ## do it by default only for bp?
                      ## tests at the beginning to see if the whole thing is feasible?
) {
  ## Checks:
  match.arg(de_type, choices = c("up_and_down", "up", "down"), several.ok = FALSE)

  # Check if there is any input at all
  if (is.null(c(de_genes, bg_genes, dds, res_de))) {
    stop(
      "Please provide one of the following forms of input: \n",
      "A vector of differentially expressed genes and a vector of backgoud genes. \n",
      "A DEseqDataset  and a DESeqResults generated by DESeq2."
    )
  }

  # Check if there only a res_de is given
  if (!is.null(res_de) & is.null(dds)) {
    stop("Please also provide a DESeq2Dataset (dds) object.")
  }

  # check if only dds is given
  if (!is.null(dds) & is.null(res_de)) {
    stop("Please also provide a DESeq2 result object.")
  }

  # Check if there only a bg but no de genes is given
  if (!is.null(bg_genes) & is.null(de_genes)) {
    stop("Please also provide a vector of  differentially expressed genes.")
  }

  # Check if there is only a de but no bg genes given
  if (!is.null(de_genes) & is.null(bg_genes)) {
    stop("Please also provide  a vector of background genes.")
  }

  if ((de_type == "up" | de_type == "down") && !is.null(de_genes)) {
    stop(
      "The argument de_type can only be used if a dds and a res_de object are provided:\n",
      "please either provide these objects or if you want to work with gene vectors set de_type to: 'up_and_down'"
    )
  }

  if (!is.null(res_de) && !is.null(dds)) {
    # Check if the inputs are the correct type

    if (!is(dds, "DESeqDataSet")) {
      stop("The provided `dds` is not a DESeqDataSet object, please check your input parameters.")
    }

    if (!is(res_de, "DESeqResults")) {
      stop("The provided `res_de` is not a DESeqResults object, please check your input parameters.")
    }

    # checking that results and dds are related
    ## at least a subset of dds should be in res
    if (!all(rownames(res_de) %in% rownames(dds))) {
      warning(
        "It is likely that the provided `dds` and `res_de` objects are not related ",
        "to the same dataset (the row names of the results are not all in the dds). ",
        "Are you sure you want to proceed?"
      )
    }

    # Check if DESeq was run on the dds
    if (!"results" %in% mcols(mcols(dds))$type) {
      stop("I couldn't find results in your dds. You should first run DESeq2::DESeq() on your dds.")
    }

    if (de_type == "up_and_down") {
      res_de_subset <- deseqresult2df(res_de, FDR = 0.05)
    } else if (de_type == "up") {
      res_de_subset <- deseqresult2df(res_de, FDR = 0.05)
      res_de_subset <- res_de_subset[res_de_subset$log2FoldChange >= 0, ]
    } else if (de_type == "down") {
      res_de_subset <- deseqresult2df(res_de, FDR = 0.05)
      res_de_subset <- res_de_subset[res_de_subset$log2FoldChange <= 0, ]
    }

    # in example top 100 but this makes more sense no?
    de_genes <- res_de_subset$id

    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]
    }
    bg_genes <- rownames(dds)[rowSums(counts(dds)) > min_counts]
    if (verbose) {
      message(
        "Your dataset has ",
        nrow(res_de_subset),
        " DE genes. You selected ",
        length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / nrow(res_de_subset)) * 100), # sprintf format with 2 decimal places
        ") genes. You analysed all ",
        de_type,
        "-regulated genes"
      )
    }
  } else if (!is.null(c(bg_genes, de_genes))) {
    all_de <- length(de_genes)

    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]
    }
    if (verbose) {
      message(
        "Your dataset has ",
        all_de,
        " DE genes.You selected ",
        length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / all_de) * 100), # sprintf format with 2 decimal places
        ") genes. You analysed all ",
        de_type,
        "-regulated genes"
      )
    }
  }

  # creating the vectors
  gene.vector <- as.integer(bg_genes %in% de_genes)
  names(gene.vector) <- bg_genes

  fdr <- FDR_GO_cutoff

  pwf <- nullp(DEgenes = gene.vector, genome = genome, id = id, plot.fit = FALSE)

  goseq_out <- goseq(pwf, genome = genome, id = id, test.cats = testCats)

  goseq_out$p.adj <- p.adjust(goseq_out$over_represented_pvalue, method = "BH")

  # removing nTop to allow for more evaluated genes
  # to reduce the load for adding the genes
  # goseq_out <- goseq_out[seq_len(nTop), ]

  if (add_gene_to_terms) {
    # for adding the gene ids/names...
    gene2cat <- getgo(de_genes, genome = genome, id = id, fetch.cats = testCats)
    names(gene2cat) <- de_genes

    reversemap <- function(map) # as in goseq
    {
      tmp <- unlist(map, use.names = FALSE)
      names(tmp) <- rep(names(map), times = as.numeric(summary(map)[, 1]))
      return(split(names(tmp), as.vector(tmp)))
    }

    cat2gene <- reversemap(gene2cat)
    # one list per GO term
    goseq_out$genes <- lapply(goseq_out$category, function(x) cat2gene[[x]])

    all_ens_ids <- unique(
      unique(unlist(goseq_out$genes))
    )

    all_genesymbols <- mapIds(get(mapping),
                              keys = all_ens_ids,
                              keytype = "ENSEMBL",
                              column = "SYMBOL",
                              multiVals = "first"
    )

    # building the lookup table
    lut_genes <- data.frame(
      gene_id = all_ens_ids,
      gene_name = all_genesymbols,
      row.names = all_ens_ids
    )

    # and also TODO: do this only if genes are not already symbols
    goseq_out$genesymbols <- lapply(goseq_out$genes, function(x) {
      sort(lut_genes[x, "gene_name"])
    })

    # coerce to char
    goseq_out$genes <- unlist(lapply(goseq_out$genes, function(arg) paste(arg, collapse = ",")))
    # coerce to char
    goseq_out$genesymbols <- unlist(lapply(goseq_out$genesymbols, function(arg) paste(arg, collapse = ",")))
  }

  return(goseq_out)
}



#' Extract functional terms enriched in the DE genes, based on clusterProfiler
#'
#' A wrapper for extracting functional GO terms enriched in a list of (DE) genes,
#' based on the algorithm and the implementation in the clusterProfiler package
#'
#' Note: the feature length retrieval is based on the \code{\link{enrichGO}} function
#'
#' @param res_de A DESeqResults object created using \code{DESeq2}
#' @param dds A DESeqDataset object created using \code{DESeq2}
#' @param de_genes A vector of (differentially expressed) genes
#' @param bg_genes A vector of background genes, e.g. all (expressed) genes in the assays
#' @param top_de numeric, how many of the top differentially expressed genes to use for the enrichment analysis.
#'  Attempts to reduce redundancy. Assumes the data is sorted by padj (default in DESeq2).
#' @param min_counts numeric, min number of counts a gene needs to have to be included
#' in the geneset that the de genes are compared to. Default is 0, recommended only for advanced users.
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param de_type One of: 'up', 'down', or 'up_and_down' Which genes to use for GOterm calculations:
#' @param keyType Gene format to input into enrichGO from clusterProfiler. If res_de and dds are used use "SYMBOL" for more
#' information check the enrichGO documentation
#' @param verbose Logical, whether to add messages telling the user which steps were taken
#' @param ... Further parameters to use for the go_enrich function from \code{clusterProfiler}
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#' @export
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom clusterProfiler enrichGO
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @importFrom DESeq2 counts
#'
#' @family Enrichment functions
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' data(gse, package = "macrophage")
#'
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#' dds_macrophage <- DESeq(dds_macrophage)
#' data(res_de_macrophage, package = "mosdef")
#'
#' library("AnnotationDbi")
#' library("org.Hs.eg.db")
#' library("clusterProfiler")
#' CluProde_macrophage <- cluproTable(
#'   res_de = res_macrophage_IFNg_vs_naive,
#'   dds = dds_macrophage,
#'   mapping = "org.Hs.eg.db"
#' )
run_cluPro <- function(res_de = NULL,
                       dds = NULL,
                       de_genes = NULL,
                       bg_genes = NULL,
                       top_de = NULL,
                       min_counts = 0,
                       mapping = "org.Hs.eg.db",
                       de_type = "up_and_down",
                       keyType = "SYMBOL",
                       verbose = TRUE,
                       ...) {
  if (!is.null(res_de) & !is.null(dds)) {
    keyType <- "SYMBOL"
    # Making sure that if res_de and dds are provided the keytype can't be overwritten to be ENTREZ
    # or something similar as the genevectors this functions creates from the res_de and dds use SYMBOLS
    # therefore enrichGO would not run with a keytype that was set wrong before
  }

  ## Checks:

  # Check if de_type is correct
  # if(!(de_type %in% c("up_and_down","up", "down")))

  match.arg(de_type, choices = c("up_and_down", "up", "down"), several.ok = FALSE)


  # Check if there is any input at all
  if (is.null(c(de_genes, bg_genes, dds, res_de))) {
    stop(
      "Please provide one of the following forms of input: \n",
      "A vector of differentially expressed genes and a vector of backgoud genes. \n",
      "A DEseqDataset  and a DESeqResults generated by DESeq2."
    )
  }

  # Check if there only a res_de is given
  if (!is.null(res_de) & is.null(dds)) {
    stop("Please also provide a DESeq2Dataset (dds) object.")
  }

  # check if only dds is given
  if (!is.null(dds) & is.null(res_de)) {
    stop("Please also provide a DESeq2 result object.")
  }

  # Check if there only a bg but no de genes is given
  if (!is.null(bg_genes) & is.null(de_genes)) {
    stop("Please also provide a vector of  differentially expressed genes.")
  }

  # Check if there is only a de but no bg genes given
  if (!is.null(de_genes) & is.null(bg_genes)) {
    stop("Please also provide  a vector of background genes.")
  }

  # checking the additional topGO_method2

  if ((de_type == "up" | de_type == "down") && !is.null(de_genes)) {
    stop(
      "The argument de_type can only be used if a dds and a res_de object are provided:\n",
      "please either provide these objects or if you want to work with gene vectors set de_type to: 'up_and_down'"
    )
  }

  annot_to_map_to <- get(mapping)

  if (!is.null(res_de) && !is.null(dds)) {
    # Check if the inputs are the correct type

    if (!is(dds, "DESeqDataSet")) {
      stop("The provided `dds` is not a DESeqDataSet object, please check your input parameters.")
    }

    if (!is(res_de, "DESeqResults")) {
      stop("The provided `res_de` is not a DESeqResults object, please check your input parameters.")
    }

    # checking that results and dds are related
    ## at least a subset of dds should be in res
    if (!all(rownames(res_de) %in% rownames(dds))) {
      warning(
        "It is likely that the provided `dds` and `res_de` objects are not related ",
        "to the same dataset (the row names of the results are not all in the dds). ",
        "Are you sure you want to proceed?"
      )
    }

    # Check if DESeq was run on the dds
    if (!"results" %in% mcols(mcols(dds))$type) {
      stop("I couldn't find results in your dds. You should first run DESeq2::DESeq() on your dds.")
    }

    res_de$symbol <- AnnotationDbi::mapIds(annot_to_map_to,
      keys = row.names(res_de),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    res_de$entrez <- AnnotationDbi::mapIds(annot_to_map_to,
      keys = row.names(res_de),
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    resOrdered <- as.data.frame(res_de[order(res_de$padj), ])

    if (de_type == "up_and_down") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
    } else if (de_type == "up") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
      de_df <- de_df[de_df$log2FoldChange >= 0, ]
    } else if (de_type == "down") {
      de_df <- resOrdered[resOrdered$padj <= .05 & !is.na(resOrdered$padj), ]
      de_df <- de_df[de_df$log2FoldChange <= 0, ]
    }

    de_genes <- de_df$symbol
    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]
    }
    bg_ids <- rownames(dds)[rowSums(counts(dds)) > min_counts]
    bg_genes <- AnnotationDbi::mapIds(annot_to_map_to,
      keys = bg_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    if (verbose) {
      message(
        "Your dataset has ",
        nrow(de_df),
        " DE genes. You selected ",
        length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / nrow(de_df)) * 100), # sprintf format with 2 decimal places
        ") genes. You analysed all ",
        de_type,
        "-regulated genes"
      )
    }
  } else if (!is.null(c(bg_genes, de_genes))) {
    all_de <- length(de_genes)

    if (!is.null(top_de)) {
      top_de <- min(top_de, length(de_genes))
      de_genes <- de_genes[seq_len(top_de)]

      if (verbose) {
        message(
          "Your dataset has ",
          all_de,
          " DE genes.You selected ",
          length(de_genes), " (", sprintf("%.2f%%", (length(de_genes) / all_de) * 100), # sprintf format with 2 decimal places
          ") genes. You analysed all ",
          de_type,
          "-regulated genes"
        )
      }
    }
  }

  res_enrich <- enrichGO(
    gene = de_genes,
    universe = bg_genes,
    OrgDb = mapping,
    keyType = keyType,
    ...
  )

  return(res_enrich)
}
