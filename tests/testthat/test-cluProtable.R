test_that("enrich_result is created", {
   CluProde_airway <- cluproTable(
     res_de = res_airway,
     dds = dds_airway,
     mapping = "org.Hs.eg.db"
   )
  expect_s4_class(CluProde_airway, "enrichResult")
  
  
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

test_that("enrich_result is created  when using vectors not res_De/dds", {
  CluProde_airway_vec <- cluproTable(
    de_genes =  myde,
    bg_genes  = myassayed,
    mapping = "org.Hs.eg.db",
    keyType = "ENSEMBL"
  )
  expect_s4_class(CluProde_airway_vec, "enrichResult")
  
})
  
  