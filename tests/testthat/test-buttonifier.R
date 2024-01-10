test_that("Output is a list for DT", {
  test <- buttonifier(airway_df)
  expect_type(test, "list")
  
})

test_that("Output is a list", {
  test <- buttonifier(airway_df, output_format = "DF")
  expect_s3_class(test, "data.frame")
  
})

test_that("All columns are created", {
  test <- buttonifier(airway_df, new_cols = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"))
  expect_type(test, "list")
  
})

test_that("All columns are created", {
  test <- buttonifier(airway_df, new_cols = c("GC", "NCBI", "GTEX", "UNIPROT", "dbPTM", "HPA"),
                      output_format = "DF")
  expect_s3_class(test, "data.frame")
  })
test_that("Errors are triggered correctly", {
  
  expect_error(buttonifier(airway_df, col_to_use = "symbol"))
  expect_error(buttonifier(airway_df, new_cols = "Genes"))
  expect_warning(buttonifier(airway_df, new_cols = c( "GC", "NCBI", "Genes")))
  
  
  
})