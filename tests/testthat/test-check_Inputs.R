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

  
  
    
})
