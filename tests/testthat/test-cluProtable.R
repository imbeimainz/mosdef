test_that("enrich results are created and top_de works properly", {
  CluProde_macrophage_topde <- cluproTable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s4_class(CluProde_macrophage_topde, "enrichResult")
  
  CluProde_macrophage_topde_vectors <- cluproTable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400,
    keyType = "ENSEMBL"
  )
  expect_s4_class(CluProde_macrophage_topde_vectors, "enrichResult")
})


test_that("enrich_result is created only for up or down -regulated genes", {
  CluProde_macrophage_up <- cluproTable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s4_class(CluProde_macrophage_up, "enrichResult")

  CluProde_macrophage_down <- cluproTable(
    res_de = res_macrophage_IFNg_vs_naive,
    dds = dds_macrophage,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s4_class(CluProde_macrophage_down, "enrichResult")
})
