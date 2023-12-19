# When running I got this error:
#   Can't find hg38/ensGene length data in genLenDataBase...
# A TxDb annotation package exists for hg38. Consider installing it to get the gene lengths.
# Trying to download from UCSC. This might take a couple of minutes. 
# Error in value[[3L]](cond) : 
#   Length information for genome hg38 and gene ID ensGene is not available. You will have to specify bias.data manually.
# fixed by installing https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html 
#not sure if that is correct but we might need to fix that


# TODO remove Ntop parameter, problem: you get multiple thousands of this message:
# 'select()' returned 1:1 mapping between keys and columns
# and then it fails with the error:
# Error in mapIds_base(x, keys, column, keytype, ..., multiVals = multiVals) : 
# apIds must have at least one key to match against.
# I assume for one of the background genes there is no database entry which then crashes everything


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
#' @param genome A string identifying the genome that genes refer to, as in the
#' \code{\link{goseq}} function
#' @param id A string identifying the gene identifier used by genes, as in the
#' \code{\link{goseq}} function
#' @param de_type One of: 'up', 'down', or 'up_and_down' Which genes to use for GOterm calculations:
#'  upregulated, downregulated or both
#' @param testCats A vector specifying which categories to test for over representation amongst DE genes - can be any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG"
#' @param FDR_GO_cutoff Numeric value for subsetting the results
#' @param nTop Number of categories to extract, and optionally process for adding
#' genes to the respective terms
#' @param orgDbPkg Character string, named as the \code{org.XX.eg.db}
#' package which should be available in Bioconductor
#' @param add_gene_to_terms Logical, whether to add a column with all genes annotated
#' to each GO term
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#' @export
#'
#' @examples
#'
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'   colData = colData(airway),
#'   design = ~ cell + dex
#' )
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#'
#' res_subset <- deseqresult2df(res_airway)[1:100, ]
#' myde <- res_subset$id
#' myassayed <- rownames(res_airway)
#' \dontrun{
#' mygo <- goseqTable(de_genes = myde,
#'   bg_genes = myassayed,
#'   testCats = "GO:BP",
#'   add_gene_to_terms = FALSE
#' )
#' head(mygo)
#' }
#'
goseqTable <- function(res_de = NULL,
                       dds = NULL,
                       de_genes = NULL, # Differentially expressed genes
                       bg_genes = NULL, # background genes, normally = rownames(cds) or filtering to genes
                       #  with at least 1 read - could also be ls(org.Mm.egGO)
                       genome = "hg38",
                       id = "ensGene",
                       de_type = "up_and_down",
                       testCats = c("GO:BP", "GO:MF", "GO:CC"),
                       FDR_GO_cutoff = 1,
                       #nTop = 200,
                       orgDbPkg = "org.Hs.eg.db",
                       # testKegg=TRUE,
                       # keggObject=mapPathwayToName("mmu"), # need the dedicated function!!
                       # writeOutput=FALSE,
                       add_gene_to_terms = TRUE # ,
                       # outputFiles_goseq="",outputFiles_goseq_kegg=""
                       ## TODO TODO: bring back in action the function
                       ## add genes annotated to each term
                       ## do it by default only for bp?
                       ## tests at the beginning to see if the whole thing is feasible?
) {
  library(goseq) # for now (again not sure how the dependencies work)
  library(GO.db)
  library(DESeq2)
  
  ##Checks:
  
  match.arg(de_type, choices = c("up_and_down","up", "down"), several.ok = FALSE)
  
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
    
    if((de_type == "up" | de_type == "down")&& !is.null(de_genes))
      stop("The argument de_type can only be used if a dds and a res_de object are provided:\n",
           "please either provide these objects or if you want to work with gene vectors set de_type to: 'up_and_down'")
    
    
        
    if(de_type == "up_and_down"){
      
      res_de_subset <- deseqresult2df(res_de)[res_de$padj <= 0.05, ] 
      
    } else if( de_type == "up") {
      
      res_de_subset <- deseqresult2df(res_de)[res_de$padj <= 0.05, ] 
      res_de_subset <- res_de_subset[res_de_subset$log2FoldChange >= 0,]
      
    } else if( de_type == "down"){
      
      res_de_subset <- deseqresult2df(res_de)[res_de$padj <= 0.05, ] 
      res_de_subset <- res_de_subset[res_de_subset$log2FoldChange <= 0,]
      
    }
   
    
    # in example top 100 but this makes more sense no?
    de_genes <- res_de_subset$id
    bg_genes <- rownames(res_de)
    
    # creating the vectors
    gene.vector <- as.integer(bg_genes %in% de_genes)
    names(gene.vector) <- bg_genes
    
  } else if(!is.null(c(bg_genes,de_genes))){
    
    # creating the vectors
    gene.vector <- as.integer(bg_genes %in% de_genes)
    names(gene.vector) <- bg_genes
  }
  
  
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
    goseq_out$genes <- sapply(goseq_out$category, function(x) cat2gene[[x]])
    
    # TODO: replace identifiers/annotaions!!!
    ## and also TODO: do this only if genes are not already symbols
    goseq_out$genesymbols <- sapply(goseq_out$genes, function(x) sort(AnnotationDbi::mapIds(get(orgDbPkg), keys = x, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")))
    # coerce to char
    goseq_out$genes <- unlist(lapply(goseq_out$genes, function(arg) paste(arg, collapse = ",")))
    # coerce to char
    goseq_out$genesymbols <- unlist(lapply(goseq_out$genesymbols, function(arg) paste(arg, collapse = ",")))
  }
  
  return(goseq_out)
}
