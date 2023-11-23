#' Create a button for genesymbols in an RMD
#' 
#' A function to turn Gene Symbols into buttons in an RMD linking to various Portals for further info
#' about these genes. 
#' Current supported portals are: Genecards, NCBI, GTEX, Uniprot, dbPTM, Human Protein Atlas
#'
#' @param df A dataframe with at least on coloumn with gene Symbols named: SYMBOL
#' @param new_cols At least one of: "GC", "NCBI", "GTEX", "UNIPROT"; "dbPTM", "HPA" 
#' 
#'
#' @return
#' @export
#'
#' @examples

#' data("gse", package = "macrophage")
#' 
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
# changing the ids to Ensembl instead of the Gencode used in the object
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#' dds_macrophage
#' dds_macrophage <- DESeq(dds_macrophage)
#' 
#' res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
#'                                         contrast = c("condition", "IFNg", "naive"),
#'                                         lfcThreshold = 1, alpha = 0.05)
#' res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL
#' res_df <- as.data.frame(res_macrophage_IFNg_vs_naive@listData)
#' res_df <-res_df[1:100,]
#' res_df <- buttons(res_df)
#' DT::datatable(res_df)

gene_symbol_buttons <- function(df, new_cols = c("GC", "UNIPROT")){
  
  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
  
  val <- df$SYMBOL
  
  #GeneCards
  if( "GC" %in% new_cols){
    for (i in 1:length(df$SYMBOL)) {
      
      
      df$SYMBOL_GC[i] <- sprintf(
        '<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn     btn-primary" style = "%s">%s</a>',
        val[i],
        .actionbutton_biocstyle,
        val[i]
      )
    }
  }
  
  #NCBI
  if( "NCBI" %in% new_cols){
    
    
    for (i in 1:length(df$SYMBOL)) {
      
      df$SYMBOL_NCBI[i] <- paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',val[i],'[sym]" target="_blank" class="btn btn-primary">',val[i],'</a>')
    }
  }
  
  #GTEx 
  if( "GTEX" %in% new_cols){
    
    
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_GTEX[i] <- sprintf(
        '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
        val[i],
        .actionbutton_biocstyle,
        val[i]
      )
    }
  }
  
  #Uniprot
  if( "UNIPROT" %in% new_cols){
    
    
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_UniProt[i] <- sprintf('<a href = "https://www.uniprot.org/uniprot/?query=%s&sort=score" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
                                      val[i],
                                      .actionbutton_biocstyle,
                                      val[i])
    }
  }
  
  #dbPTM 
  if( "dbPTM" %in% new_cols){
    
    
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_dbPTM[i] <- sprintf('<a href = "http://dbptm.mbc.nctu.edu.tw/info.php?id=%s_HUMAN" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-edit"></i>%s</a>',
                                    val[i],
                                    .actionbutton_biocstyle,
                                    val[i])
    }
  }
  
  #Human protein atlas Portal
  if( "HPA" %in% new_cols){
    
    
    for (i in 1:length(df$SYMBOL)) {
      df$SYMBOL_HPA[i] <- sprintf('<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
                                  val[i],
                                  .actionbutton_biocstyle,
                                  val[i])
    }
  }
  
  df <- df %>%
    select(-SYMBOL)
  return(df)
}
