test_that("Enrichresult is created and top_de works properly", {
  CluProde_airway_topde <- cluproTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s4_class(CluProde_airway_topde, "enrichResult")
  
  CluProde_airway_topde_vectors <- cluproTable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400,
    keyType = "ENSEMBL"
  )
  expect_s4_class(CluProde_airway_topde_vectors, "enrichResult")
})

test_that("enrich_result is created only for up or down -regulated genes", {
  CluProde_airway_up <- cluproTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s4_class(CluProde_airway_up, "enrichResult")

  CluProde_airway_down <- cluproTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s4_class(CluProde_airway_down, "enrichResult")
})



