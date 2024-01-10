context("Input main parameter checking works")

test_that("Early fails are triggered", {
  
  #res_de is a DESeqResults
  expect_error({
    cluproTable(
      res_de = myde,
      dds = dds_airway,
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    de_volcano(res_de = myde,
               L2FC_cutoff = 1,
               labeled_genes = 20,
               mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    deseqresult2df(myde) 
  })
  
  expect_error({
    deseqresult2DEgenes(myde) 
  })
  
  expect_error({
    deseqresult2tbl(myde) 
  })
  
  expect_error({
    goseqTable(
      res_de = myde,
      dds = dds_airway,
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    plot_ma(res_de = myde, 
            FDR = 0.05,
            hlines = 1)
  })
  
  expect_error({
    topGOtable(res_de = myde,
               dds = dds_airway,
               ontology = "BP",
               mapping = "org.Hs.eg.db",
               geneID = "symbol"
    )
  })
})
test_that("check if dds are dds and vecs are vecs", {
  expect_error({
    cluproTable(
      res_de = res_airway,
      dds = myassayed,
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    gene_plot(res_airway,
              gene = "ENSG00000125347",
              intgroup = "condition",
              annotation_obj = anno_df
    )
  })
  
  expect_error({
    goseqTable(
      res_de = res_airway,
      dds = myassayed,
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    topGOtable(res_de = res_airway,
               dds = myassayed,
               ontology = "BP",
               mapping = "org.Hs.eg.db",
               geneID = "symbol"
    )
  })
  
  # check if vectors are actually vectors
  expect_error({
    cluproTable(
      de_genes = res_airway,
      bg_genes = dds_airway,
      mapping = "org.Hs.eg.db"
    )
  })
  
  
  expect_error({
    goseqTable(
      de_genes = res_airway,
      bg_genes = dds_airway,
      mapping = "org.Hs.eg.db"
    )
  })
  
  
  expect_error({
    topGOtable(de_genes = res_airway,
               bg_genes = dds_airway,
               ontology = "BP",
               mapping = "org.Hs.eg.db",
               geneID = "symbol"
    )
  })
})

test_that("Error is thrown with insufficient input", {
  #check insufficient input
  
  expect_error({
    cluproTable(
      mapping = "org.Hs.eg.db"
    )
  })
  
  
  
  expect_error({
    goseqTable(
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_error({
    topGOtable(
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})

test_that("Check if de_type is correct", {
  
  expect_error({
    cluproTable(
      res_de = res_airway,
      dds = dds_airway,
      mapping = "org.Hs.eg.db",
      de_type = "all"
    )
  })
  
  expect_error({
    goseqTable(
      res_de = res_airway,
      dds = dds_airway,
      mapping = "org.Hs.eg.db",
      de_type = "all"
    )
  })
  
  expect_error({
    topGOtable(res_de = res_airway,
               dds = dds_airway,
               ontology = "BP",
               mapping = "org.Hs.eg.db",
               geneID = "symbol",
               de_type = "all"
    )
  })
  
  
})

test_that("res_de and dds are related", {
  
  
  expect_warning({
    cluproTable(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_airway,
      mapping = "org.Hs.eg.db"
    )
  })
  
  expect_warning({
    goseqTable(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_airway,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })
  
  expect_warning({
    topGOtable(res_de = res_macrophage_IFNg_vs_naive,
               dds = dds_airway,
               ontology = "BP",
               mapping = "org.Hs.eg.db",
               geneID = "symbol"
    )
  })
  
  
})










