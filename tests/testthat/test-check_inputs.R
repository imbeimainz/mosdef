test_that("Early fails are triggered", {
  # res_de is a DESeqResults
  expect_error({
    run_cluPro(
      res_de = myde,
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    de_volcano(
      res_de = myde,
      L2FC_cutoff = 1,
      labeled_genes = 20,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    deseqresult2df(myde)
  })
  expect_error({
    signature_volcano(myde)
  })

  expect_error({
    deseqresult2DEgenes(myde)
  })

  expect_error({
    run_goseq(
      res_de = myde,
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    plot_ma(
      res_de = myde,
      FDR = 0.05,
      hlines = 1
    )
  })

  expect_error({
    run_topGO(
      res_de = myde,
      dds = dds_macrophage,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})
test_that("check if dds are dds and vecs are vecs", {
  expect_error({
    run_cluPro(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = myassayed,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    gene_plot(res_macrophage_IFNg_vs_naive,
      gene = "ENSG00000125347",
      intgroup = "condition",
      annotation_obj = anno_df
    )
  })

  expect_error({
    run_goseq(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = myassayed,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_topGO(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = myassayed,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })

  # check if vectors are actually vectors
  expect_error({
    run_cluPro(
      de_genes = res_macrophage_IFNg_vs_naive,
      bg_genes = dds_macrophage,
      mapping = "org.Hs.eg.db"
    )
  })


  expect_error({
    run_goseq(
      de_genes = res_macrophage_IFNg_vs_naive,
      bg_genes = dds_macrophage,
      mapping = "org.Hs.eg.db"
    )
  })


  expect_error({
    run_topGO(
      de_genes = res_macrophage_IFNg_vs_naive,
      bg_genes = dds_macrophage,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})

test_that("Error is thrown with insufficient input", {
  # check insufficient input

  expect_error({
    run_cluPro(
      mapping = "org.Hs.eg.db"
    )
  })



  expect_error({
    run_goseq(
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_topGO(
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})

test_that("Check if de_type is correct", {
  expect_error({
    run_cluPro(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db",
      de_type = "all"
    )
  })

  expect_error({
    run_goseq(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db",
      de_type = "all"
    )
  })

  expect_error({
    run_topGO(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol",
      de_type = "all"
    )
  })
})

test_that("res_de and dds are related", {
  expect_warning(
    {
      run_cluPro(
        res_de = res_macrophage_IFNg_vs_naive,
        dds = dds_airway,
        mapping = "org.Hs.eg.db"
      )
    },
    "not related"
  )
})

test_that("res_de and dds are related", {
  expect_warning(
    {
      run_goseq(
        res_de = res_macrophage_IFNg_vs_naive,
        dds = dds_airway,
        mapping = "org.Hs.eg.db",
        add_gene_to_terms = FALSE
      )
    },
    "not related"
  )
})

test_that("res_de and dds are related", {
  expect_warning(
    {
      run_topGO(
        res_de = res_macrophage_IFNg_vs_naive,
        dds = dds_airway,
        ontology = "BP",
        mapping = "org.Hs.eg.db",
        geneID = "symbol"
      )
    },
    "not related"
  )
})

test_that("DESeq was run on the dds", {
  expect_error({
    run_cluPro(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage_nodeseq,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_goseq(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage_nodeseq,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })

  expect_error({
    run_topGO(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage_nodeseq,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})

test_that("Errors are thrown if only one of two needed inputs is provided", {
  # res_de missing


  expect_error({
    run_cluPro(
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_goseq(
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })

  expect_error({
    run_topGO(
      dds = dds_macrophage,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })

  # dds is missing

  expect_error({
    run_cluPro(
      res_de = res_macrophage_IFNg_vs_naive,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_goseq(
      res_de = res_macrophage_IFNg_vs_naive,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })

  expect_error({
    run_topGO(
      res_de = res_macrophage_IFNg_vs_naive,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })

  # de_genes is missing


  expect_error({
    run_cluPro(
      bg_genes = myassayed,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_goseq(
      bg_genes = myassayed,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })

  expect_error({
    run_topGO(
      bg_genes = myassayed,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })

  # bg_genes is missing


  expect_error({
    run_cluPro(
      de_genes = myde,
      mapping = "org.Hs.eg.db"
    )
  })

  expect_error({
    run_goseq(
      de_genes = myde,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE
    )
  })

  expect_error({
    run_topGO(
      de_genes = myde,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol"
    )
  })
})

test_that("de_type can not be used with vectors to avoid confusion", {
  expect_error({
    run_cluPro(
      de_genes = myde,
      bg_genes = myassayed,
      mapping = "org.Hs.eg.db",
      de_type = "up"
    )
  })

  expect_error({
    run_goseq(
      de_genes = myde,
      bg_genes = myassayed,
      mapping = "org.Hs.eg.db",
      add_gene_to_terms = FALSE,
      de_type = "up"
    )
  })

  expect_error({
    run_topGO(
      de_genes = myde,
      bg_genes = myassayed,
      ontology = "BP",
      mapping = "org.Hs.eg.db",
      geneID = "symbol",
      de_type = "up"
    )
  })
})
