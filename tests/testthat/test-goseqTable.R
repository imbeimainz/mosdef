test_that("enrich_result is created", {
  goseqde_airway <- goseqTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db"
  )
  expect_s3_class(goseqde_airway, "data.frame")
})

test_that("enrich_result is created only for up or down -regulated genes", {
  goseqde_airway_up <- goseqTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s3_class(goseqde_airway_up, "data.frame")

  goseqde_airway_down <- goseqTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s3_class(goseqde_airway_down, "data.frame")
})

test_that("enrich_result is created  when using vectors not res_De/dds", {
  goseqde_airway_vec <- goseqTable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db"
  )
  expect_s3_class(goseqde_airway_vec, "data.frame")
})

test_that("top_de works propperly", {
  goseqde_airway_topde <- goseqTable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(goseqde_airway_topde, "data.frame")

  goseqde_airway_topde_vectors <- goseqTable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(goseqde_airway_topde_vectors, "data.frame")
})
