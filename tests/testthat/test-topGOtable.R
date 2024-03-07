test_that("enrich_result is created and top_de works", {
  res_enrich_macrophage_topGO_topde <- topGOtable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(res_enrich_macrophage_topGO_topde, "data.frame")
  
  res_enrich_macrophage_topGO_topde_vectors <- topGOtable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400,
    geneID = "ensembl"
  )
  expect_s3_class(res_enrich_macrophage_topGO_topde_vectors, "data.frame")
})

test_that("enrich_result is created only for up or down -regulated genes", {
  res_enrich_macrophage_topGO_up <- topGOtable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s3_class(res_enrich_macrophage_topGO_up, "data.frame")

  res_enrich_macrophage_topGO_down <- topGOtable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s3_class(res_enrich_macrophage_topGO_down, "data.frame")
})


test_that("Error is thrown if wrong topGO method is used", {
  expect_error(
    topGOtable(
      res_de = res_macrophage_IFNg_vs_naive,
      dds = dds_macrophage,
      mapping = "org.Hs.eg.db",
      topGO_method2 = "test"
    )
  )
})
test_that("do_paj works", {
  res_enrich_macrophage_topGO_padj <- topGOtable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    do_padj = TRUE
  )
  expect_s3_class(res_enrich_macrophage_topGO_padj, "data.frame")
})
