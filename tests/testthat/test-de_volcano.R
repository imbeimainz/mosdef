test_that(" ggplot object is created", {
  p <- de_volcano(res_airway,
    mapping = "org.Hs.eg.db"
  )
  expect_s3_class(p, "gg")

  # LFC Cutogg and labeled genes can be used
  p_L2FC <- de_volcano(res_airway,
    mapping = "org.Hs.eg.db",
    L2FC_cutoff = 2
  )
  expect_s3_class(p_L2FC, "gg")

  p_label <- de_volcano(res_airway,
    mapping = "org.Hs.eg.db",
    labeled_genes = 30
  )
  expect_s3_class(p_label, "gg")
})
