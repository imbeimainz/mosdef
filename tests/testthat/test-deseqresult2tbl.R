test_that("Result df is created", {
  result_tbl <- deseqresult2tbl(res_macrophage_IFNg_vs_naive)
  expect_s3_class(result_tbl, "data.frame")
})
