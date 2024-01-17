test_that("ggplot is created", {
  
  p <- plot_ma(res_airway, FDR = 0.05, hlines = 1)
  expect_s3_class(p, "gg")
  
})

test_that("intgenes can be used",code = {
  
  p_intgenes <- plot_ma(res_airway,
    FDR = 0.1,
    intgenes = c(
      "ENSG00000103196", # CRISPLD2
      "ENSG00000120129", # DUSP1
      "ENSG00000163884", # KLF15
      "ENSG00000179094" # PER1
    )
  )
  expect_s3_class(p_intgenes, "gg")
  
  
})

test_that("Other parameters can be used if not at default value", {
  
  p_other <- plot_ma(res_airway,
                     point_alpha = 0.02,
                     draw_y0 = FALSE,
                     hlines = 2,
                     title = "test",
                     xlab = "mean of normalized counts - log10 scale",
                     ylim = 5,
                     add_rug = FALSE,
                     intgenes_color = "steelblue",
                     labels_intgenes = FALSE,
                     labels_repel = FALSE)
  
  expect_s3_class(p_other, "gg")
  
   
})
