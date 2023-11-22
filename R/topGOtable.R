
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
topGOtable <- function(res_de = NULL,                  # Differentially expressed genes
                       dds = NULL,                 # background genes
                       de_genes = NULL,
                       bg_genes = NULL,
                       ontology = "BP",            # could use also "MF"
                       annot = annFUN.org,       # parameters for creating topGO object
                       mapping = "org.Mm.eg.db",
                       geneID = "symbol",       # could also beID = "entrez"
                       #topTablerows = 200,     was used to limit table size. Now just use nrows
                       fullNamesInRows = TRUE,
                       addGeneToTerms = TRUE,
                       plotGraph = FALSE, 
                       plotNodes = 10,
                       writeOutput = FALSE, 
                       outputFile = "",
                       topGO_method2 = "elim",
                       do_padj = FALSE) {
  ##Checks:
  
  # Check if there is any input at all
  if(is.null(c(de_genes,bg_genes,dds, res_de)))
    
    stop("Please provide one of the following forms of input: \n",
         "A vector of differentially expressed genes and a vector of backgoud genes. \n",
         "A DEseqDataset  and a DESeqResults generated by DESeq2.")
  
  # Check if there only a res_de is given
  if(!is.null(res_de)& is.null(dds))
    
    stop("Please also provide a DESeq2Dataset (dds) object.")
  
  #check if only dds is given
  if(!is.null(dds)& is.null(res_de))
    
    stop("Please also provide a DESeq2 result object.")
  
  # Check if there only a bg but no de genes is given
  if(!is.null(bg_genes)& is.null(de_genes))
    
    stop("Please also provide a vector of  differentially expressed genes.")
  
  # Check if there is only a de but no bg genes given
  if (!is.null(de_genes)& is.null(bg_genes))
    
    stop("Please also provide  a vector of background genes.")
  
  
  # checking the additional topGO_method2
  topgo_methods <- c("elim", "weight", "weight01", "lea", "parentchild")
  if (!(topGO_method2 %in% topgo_methods))
    stop("Please provide one of the following topGO methods in addition to the classic method:\n",
         paste(topgo_methods, collapse = ", "))
  
  #Dependencies
  library("AnnotationDbi")# for the dependencies I don't know how to set them :D
  library("topGO") # see above
  annot_to_map_to <- get(mapping)
  
  if(!is.null(res_de)&& !is.null(dds)) {
    
    #Check if the inputs are the correct type
    
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
    if(!"results" %in% mcols(mcols(dds))$type){
      stop("I couldn't find results in your dds. You should first run DESeq2::DESeq() on your dds.")
      
    }
    
    
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
    de_genes_input <- factor(as.integer(bg_symbols %in% de_symbols))
    names(de_genes_input) <- bg_symbols
    
  }   else if(!is.null(c(de_genes,bg_genes))){
    # creating the vectors
    de_genes_input <- factor(as.integer(bg_genes %in% de_genes))
    names(de_genes_input) <- bg_genes
    
  }
  
  
  # instantiating the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = de_genes_input,
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
  ### topTablerows <- min(nrow(sTab), topTablerows)
  topTablerows <- nrow(sTab)
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
