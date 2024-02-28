test_that("Result df is created", {
  result_df <- deseqresult2df(res_macrophage_IFNg_vs_naive)
  expect_s3_class(result_df, "data.frame")
})

test_that("FDR can be used instead of further subsetting later", {
  result_df <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.01)
  expect_s3_class(result_df, "data.frame")
})

test_that("Result df is created", {
  result_df <- deseqresult2DEgenes(res_macrophage_IFNg_vs_naive)
  expect_s3_class(result_df, "data.frame")
})

test_that("FDR can be used instead of further subsetting later", {
  result_df <- deseqresult2DEgenes(res_macrophage_IFNg_vs_naive, FDR = 0.01)
  expect_s3_class(result_df, "data.frame")
})
