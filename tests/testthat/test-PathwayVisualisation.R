test_that("VisualisePathways returns an informative plot for an empty min_n result", {
  pathway_df <- data.frame(
    analytes_in_pathways = c(1L, 2L),
    p_val = c(0.01, 0.02)
  )

  plot <- VisualisePathways(
    SpaMTP = NULL,
    pathway_df = pathway_df,
    min_n = 3,
    verbose = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "No pathways to visualise")
})

test_that("VisualisePathways returns an informative plot after p-value filtering", {
  pathway_df <- data.frame(
    pathway_name = "Example pathway",
    pathway_id = "example:1",
    p_val = 0.2,
    analytes_in_pathways = 3L,
    total_in_pathways = 10L,
    adduct_info = ""
  )

  plot <- VisualisePathways(
    SpaMTP = NULL,
    pathway_df = pathway_df,
    min_n = 3,
    p_val_threshold = 0.1,
    verbose = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "No pathways to visualise")
})
