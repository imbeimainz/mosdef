test_that("Result df is created", {
  result_tbl <- deseqresult2tbl(res_airway)
  expect_s3_class(result_tbl, "data.frame")
})
