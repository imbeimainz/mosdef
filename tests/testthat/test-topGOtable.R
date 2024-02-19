test_that("Enrich result is created and top_de works propperly", {
  topGOde_airway_topde <- topGOtable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    top_de = 400
  )
  expect_s3_class(topGOde_airway_topde, "data.frame")
  
  topGOde_airway_topde_vectors <- topGOtable(
    de_genes = myde,
    bg_genes = myassayed,
    mapping = "org.Hs.eg.db",
    top_de = 400,
    geneID = "ensembl"
  )
  expect_s3_class(topGOde_airway_topde_vectors, "data.frame")
})

test_that("enrich_result is created only for up or down -regulated genes", {
  topGOde_airway_up <- topGOtable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "up"
  )
  expect_s3_class(topGOde_airway_up, "data.frame")

  topGOde_airway_down <- topGOtable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    de_type = "down"
  )
  expect_s3_class(topGOde_airway_down, "data.frame")
})

test_that("Error is thrown if wrong topGO method is used", {
  expect_error(
    topGOtable(
      res_de = res_airway,
      dds = dds_airway,
      mapping = "org.Hs.eg.db",
      topGO_method2 = "test"
    )
  )
})
test_that("do_paj works", {
  topGOde_airway_padj <- topGOtable(
    res_de = res_airway,
    dds = dds_airway,
    mapping = "org.Hs.eg.db",
    do_padj = TRUE
  )
  expect_s3_class(topGOde_airway_padj, "data.frame")
})
