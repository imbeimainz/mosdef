
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
#' @param DEgenes A vector of (differentially expressed) genes
#' @param BGgenes A vector of background genes, e.g. all (expressed) genes in the assays
#' @param dds A DESeqDataset object created using \code{DESeq2} 
#' @param res_de A DESeqResults object created using \code{DESeq2} 
#' @param ontology Which Gene Ontology domain to analyze: \code{BP} (Biological Process), \code{MF} (Molecular Function), or \code{CC} (Cellular Component)
#' @param annot Which function to use for annotating genes to GO terms. Defaults to \code{annFUN.org}
#' @param annot_to_map_to A temporary parameter allowing the use of AnnotationDBs mapIds function
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param geneID Which format the genes are provided. Defaults to \code{symbol}, could also be
#' \code{entrez} or \code{ENSEMBL}
#' @param topTablerows How many rows to report before any filtering
#' @param fullNamesInRows Logical, whether to display or not the full names for the GO terms
#' @param addGeneToTerms Logical, whether to add a column with all genes annotated to each GO term
#' @param plotGraph Logical, if TRUE additionally plots a graph on the identified GO terms
#' @param plotNodes Number of nodes to plot
#' @param writeOutput Logical, if TRUE additionally writes out the result to a file
#' @param outputFile Name of the file the result should be written into
#' @param topGO_method2 Character, specifying which of the methods implemented by \code{topGO} should be used, in addition to the \code{classic} algorithm. Defaults to \code{elim}
#' @param do_padj Logical, whether to perform the adjustment on the p-values from the specific
#' topGO method, based on the FDR correction. Defaults to FALSE, since the assumption of 
#' independent hypotheses is somewhat violated by the intrinsic DAG-structure of the Gene
#' Ontology Terms
#'
#' @import topGO
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#'
#' @examples
#' library(airway)
#' library(DESeq2)
#' data(airway)
#' airway
#' dds_airway <- DESeqDataSet(airway, design= ~ cell + dex)
#' # Example, performing extraction of enriched functional categories in
#' # detected significantly expressed genes
#' \dontrun{
#' dds_airway <- DESeq(dds_airway)
#' res_airway <- results(dds_airway)
#' library("AnnotationDbi")
#' library("org.Hs.eg.db")
#' res_airway$symbol <- mapIds(org.Hs.eg.db,
#'                             keys = row.names(res_airway),
#'                             column = "SYMBOL",
#'                             keytype = "ENSEMBL",
#'                             multiVals = "first")
#' res_airway$entrez <- mapIds(org.Hs.eg.db,
#'                             keys = row.names(res_airway),
#'                             column = "ENTREZID",
#'                             keytype = "ENSEMBL",
#'                             multiVals = "first")
#' resOrdered <- as.data.frame(res_airway[order(res_airway$padj),])
#' de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
#' de_symbols <- de_df$symbol
#' bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
#' bg_symbols <- mapIds(org.Hs.eg.db,
#'                      keys = bg_ids,
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL",
#'                      multiVals = "first")
#' library(topGO)
#' topgoDE_airway <- topGOtable(de_symbols, bg_symbols,
#'                              ontology = "BP",
#'                              mapping = "org.Hs.eg.db",
#'                              geneID = "symbol")
#' }
#'
#' @export
topGOtable <- function(DEgenes = NULL,                  # Differentially expressed genes
                       BGgenes = NULL,                 # background genes
                       dds = NULL,
                       res_de = NULL,
                       ontology = "BP",            # could use also "MF"
                       annot = annFUN.org,       # parameters for creating topGO object
                       mapping = "org.Mm.eg.db",
                       annot_to_map_to = org.Mm.eg.db, # can't figure out how to turn the string from mapping 
                                                  # into the annotation needed for mapIds 
                                                  # also the name is shit
                       geneID = "symbol",       # could also beID = "entrez"
                       topTablerows = 200,
                       fullNamesInRows = TRUE,
                       addGeneToTerms = TRUE,
                       plotGraph = FALSE, 
                       plotNodes = 10,
                       writeOutput = FALSE, 
                       outputFile = "",
                       topGO_method2 = "elim",
                       do_padj = FALSE) {
  # Check if there is any input at all
  if(is.null(c(DEgenes,BGgenes,dds, res_de)))
    # and a res_de? In theory we can generate the res_de insinde topGOtable
    stop("Please provide one of the following forms of input: \n",
         "A vector of differentially expressed genes and a vector of backgoud genes \n",
         "A DEseqDataset generated by DESeq2")
  
  # Check if there only a res_de is given
  if(!is.null(res_de)& is.null(dds))
    
    stop("Please also provide a DESeq2 (dds) object")
  
  # Check if there only a bg or de genes is given
  if((!is.null(BGgenes)& is.null(DEgenes))|| (!is.null(DEgenes)& is.null(BGgenes)))
    
    stop("Please also provide both a vector of background genes and of differentially expressed genes")
  
  
  # checking the additional topGO_method2
  topgo_methods <- c("elim", "weight", "weight01", "lea", "parentchild")
  if (!(topGO_method2 %in% topgo_methods))
    stop("Please provide one of the following topGO methods in addition to the classic method:\n",
         paste(topgo_methods, collapse = ", "))
  
  #Dependencies
  library("AnnotationDbi")# for the dependencies I don't know how to set them :D
  library("topGO") # see above
  #library(package=mapping) This doesn't work but again: do we just put these packages as dependencies?
  
  
   if(!is.null(res_de)&& !is.null(dds)) {
    
    
    res_de$symbol <- AnnotationDbi::mapIds(annot_to_map_to,
                                        keys = row.names(res_de),
                                        column = "SYMBOL",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")
    res_de$entrez <- AnnotationDbi::mapIds(annot_to_map_to, 
                                        keys = row.names(res_de),
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")
    resOrdered <- as.data.frame(res_de[order(res_de$padj),])
    de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
    de_symbols <- de_df$symbol
    bg_ids <- rownames(dds)[rowSums(counts(dds)) > 0]
    bg_symbols <- mapIds(org.Hs.eg.db,
                         keys = bg_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    
    # creating the vectors
    DEgenes_input <- factor(as.integer(bg_symbols %in% de_symbols))
    names(DEgenes_input) <- bg_symbols
   } 
  
  else if (! is.null(dds)){


    dds <- DESeq(dds)
    res <- results(dds)

    res$symbol <- AnnotationDbi::mapIds(annot_to_map_to,
                                 keys = row.names(res),
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
    res$entrez <- AnnotationDbi::mapIds(annot_to_map_to, 
                                 keys = row.names(res),
                                 column = "ENTREZID",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
    resOrdered <- as.data.frame(res[order(res$padj),])
    de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
    de_symbols <- de_df$symbol
    bg_ids <- rownames(dds)[rowSums(counts(dds)) > 0]
    bg_symbols <- mapIds(org.Hs.eg.db,
                         keys = bg_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    
    # creating the vectors
    DEgenes_input <- factor(as.integer(bg_symbols %in% de_symbols))
    names(DEgenes_input) <- bg_symbols
  }

  else if(!is.null(c(DEgenes,BGgenes))){
    # creating the vectors
    DEgenes_input <- factor(as.integer(BGgenes %in% DEgenes))
    names(DEgenes_input) <- BGgenes
  }
    
 
  # instantiating the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = DEgenes_input,
                nodeSize = 10,
                annot = annot,
                mapping = mapping,
                ID = geneID)
  # performing the test(s)
  result_method2 <- runTest(GOdata, algorithm = topGO_method2, statistic = "fisher")
  resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sTab <- GenTable(GOdata,
                   p.value_method2 = result_method2,
                   p.value_classic = resultClassic,
                   orderBy = "p.value_method2",
                   ranksOf = "p.value_classic",
                   topNodes = length(score(resultClassic)))
  
  names(sTab)[which(names(sTab) == "p.value_method2")] <- paste0("p.value_", topGO_method2)
  
  sTab[["p.value_classic"]] <- as.numeric(sTab[["p.value_classic"]])
  sTab[[paste0("p.value_", topGO_method2)]] <- as.numeric(sTab[[paste0("p.value_", topGO_method2)]])
  
  # if FDR, then apply it here
  if (do_padj)
    sTab[[paste0("padj_BY_", topGO_method2)]] <- p.adjust(sTab[[paste0("p.value_", topGO_method2)]], method = "BY")
  
  # subset to specified number of rows
  topTablerows <- min(nrow(sTab), topTablerows)
  sTab <- sTab[seq_len(topTablerows), ]
  
  if (fullNamesInRows) {
    sTab$Term <- sapply(sTab$GO.ID, function(go) {
      Term(GOTERM[[go]])
    })
  }
  
  if (addGeneToTerms) {
    # adapted from an elegant one liner found here: https://support.bioconductor.org/p/65856/
    SignificantGenes <- sigGenes(GOdata)
    sTab$genes <- sapply(sTab$GO.ID, function(x) {
      genes <- genesInTerm(GOdata, x)
      tmp <- genes[[1]][genes[[1]] %in% SignificantGenes]
    })
    # coerce the list to a comma separated vector
    sTab$genes <- unlist(lapply(sTab$genes, function(arg) paste(arg, collapse = ",")))
  }
  
  # write all entries of the table
  if (writeOutput) write.table(sTab, file = outputFile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  if (plotGraph) showSigOfNodes(GOdata, topGO::score(result_method2), firstSigNodes = plotNodes, useInfo = "all")
  #   if(outputToLatex) sTabSig <- xtable(apply(sTabSig[1:15,], 2, as.character)) # take a smaller subset
  
  # and returns the significant ones # or all, like here
  return(sTab)
}
