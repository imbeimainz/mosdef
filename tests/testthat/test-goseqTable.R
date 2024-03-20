test_that("enrich results are created and top_de works properly", {
  goseqde_macrophage_topde <- run_goseq(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(goseqde_macrophage_topde, "data.frame")

  goseqde_macrophage_topde_vectors <- run_goseq(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(goseqde_macrophage_topde_vectors, "data.frame")
})


test_that("enrich_result is created only for up or down -regulated genes", {
  goseqde_macrophage_up <- run_goseq(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s3_class(goseqde_macrophage_up, "data.frame")

  goseqde_macrophage_down <- run_goseq(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s3_class(goseqde_macrophage_down, "data.frame")
})
