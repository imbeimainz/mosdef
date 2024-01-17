test_that("Result df is created", {
  
  result_df <- deseqresult2df(res_airway)
  expect_s3_class(result_df, "data.frame")
  
  
})

test_that("FDR can be used instead of further subsetting later", {
  
  result_df <- deseqresult2df(res_airway, FDR = 0.01)
  expect_s3_class(result_df, "data.frame")
  
  
})

