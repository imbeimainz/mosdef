test_that(" ggplot object is created", {
  p <- signature_volcano(res_macrophage_IFNg_vs_naive)
  expect_s3_class(p, "gg")
})

test_that("Other parameters work", {
  p_param <- signature_volcano(res_macrophage_IFNg_vs_naive, # res_enrich,
    annotation_obj = NULL,
    genelist = c(
      "ENSG00000108702",
      "ENSG00000181374",
      "ENSG00000276409"
    ),
    FDR = 0.05,
    color = "#1a81c2",
    volcano_labels = 30,
    plot_title = "Test"
  )
  expect_s3_class(p_param, "gg")
})
