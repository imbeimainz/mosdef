test_that("Functions to create buttons work", {
  goseqde_airway$button <- create_link_GO(goseqde_airway$category)
  expect_type(goseqde_airway$button, "character")
  value <- "http"
  chars <- as.character(goseqde_airway$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))

  res_subset$button <- create_link_pubmed(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))

  res_subset$button_ens <- create_link_ENS(res_subset$id)
  expect_type(res_subset$button_ens, "character")
  value <- "http"
  chars <- as.character(res_subset$button_ens[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))

  res_subset$button <- create_link_NCBI(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))

  res_subset$button <- create_link_GTEX(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))

  res_subset$button <- create_link_UniProt(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))


  res_subset$button <- create_link_dbPTM(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))



  res_subset$button <- create_link_HPA(res_subset$SYMBOL)
  expect_type(res_subset$button, "character")
  value <- "http"
  chars <- as.character(res_subset$button[1])

  # expect that http is in the strings to check if they are created properlly
  expect_true(grepl(value, chars, fixed = TRUE))
})


test_that("go_2_html works", {
  test_go <- go_2_html("GO:0009653")
  expect_s3_class(test_go, "html")

  test_go_res_enrich <- go_2_html("GO:0009653", goseqde_airway)
  expect_s3_class(test_go, "html")


  test_wrongGO <- go_2_html("yaddayadda")
  expect_s3_class(test_wrongGO, "html")
  value <- "not found"
  chars <- as.character(test_wrongGO)
  expect_true(grepl(value, chars, fixed = TRUE))
})

test_that("geneinfo_2_html works", {
  test_gene_info <- geneinfo_2_html("SPARCL1")
  expect_s3_class(test_gene_info, "html")

  test_gene_info_res <- geneinfo_2_html("SPARCL1", res_airway)
  expect_s3_class(test_gene_info_res, "html")

  expect_message(geneinfo_2_html("P53", res_airway))
})
