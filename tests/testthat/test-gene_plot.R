context("Testing function for plotting the gene expression levels")

test_that("Basic gene plot is generated", {
  p <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    annotation_obj = anno_df,
    transform = TRUE,
    labels_repel = TRUE
  )
  expect_s3_class(p, "gg")
  
  p2_noanno_normallabels_untransformed <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    transform = FALSE,
    labels_repel = FALSE
  )
  expect_s3_class(p2_noanno_normallabels_untransformed, "gg")
  
  expect_error({
    gene_plot(
      dds = dds_macrophage,
      gene = "ENSG00000285982",
      assay = "counts",
      intgroup = "factor_not_there",
      annotation_obj = anno_df,
      transform = TRUE,
      labels_repe
    )
  })
})

test_that("Enforcing plot types", {
  p_jitter <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    plot_type = "jitteronly"
  )
  p_boxplot <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    plot_type = "boxplot"
  )
  p_violin <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    plot_type = "violin"
  )
  p_sina <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    plot_type = "sina"
  )
  expect_s3_class(p_jitter, "gg")
  expect_s3_class(p_boxplot, "gg")
  expect_s3_class(p_violin, "gg")
  expect_s3_class(p_sina, "gg")
})

test_that("Data instead of plot is returned", {
  df_jitter <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    return_data = TRUE
  )
  expect_s3_class(df_jitter, "data.frame")
})

test_that("Assays are correctly accessed", {
  p_non_norm_counts <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "counts",
    intgroup = "condition",
    normalized = FALSE
  )
  expect_s3_class(p_non_norm_counts, "gg")
  p_tpm <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "abundance",
    intgroup = "condition",
    normalized = FALSE
  )
  expect_s3_class(p_tpm, "gg")
  
  p_other_assay <- gene_plot(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    assay = "avgTxLength",
    intgroup = "condition",
    normalized = FALSE
  )
  expect_s3_class(p_other_assay, "gg")
})

test_that("Extraction of expression values works", {
  df_simple <- get_expr_values(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    intgroup = "condition",
    assay = "counts"
  )
  expect_s3_class(df_simple, "data.frame")
  

  expect_error(get_expr_values(
    dds = dds_macrophage,
    gene = "ENSG00000285982",
    intgroup = "condition",
    assay = "count"
  ))
  df_unnormalized <- get_expr_values(
    dds = dds_unnormalized,
    gene = "ENSG00000285982",
    intgroup = "condition",
    assay = "counts"
  )
  expect_s3_class(df_unnormalized, "data.frame")
})
