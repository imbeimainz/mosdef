test_that("map2color works", {
  mypal <- rev(scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4
  ))
  my_vals <- res_macrophage_IFNg_vs_naive$log2FoldChange[1:20]
  m2c <- map2color(x = my_vals, pal = mypal, limits = c(-4, 4))
  m2c_nolimits <- map2color(x = my_vals, pal = mypal)
  # plot(1:20, col = m2c, pch = 20, cex = 5)
  
  expect_length(m2c, 20)
  expect_length(m2c_nolimits, 20)
})