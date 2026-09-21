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

spatial_pathway_test_object <- function() {
  counts <- rbind(
    varying = c(1, 3, 2, 5),
    constant = c(2, 2, 2, 2)
  )
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(
    counts = SeuratObject::CreateAssayObject(counts = counts), assay = "RNA"
  )
  object <- SeuratObject::SetAssayData(
    object, assay = "RNA", layer = "scale.data", new.data = counts
  )
  coords <- data.frame(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1),
                       row.names = colnames(counts))
  object[["slice1"]] <- SeuratObject::CreateFOV(
    coords, type = "centroids", assay = "RNA"
  )
  object
}

test_that("spatial pathways without measured features return a labelled panel", {
  object <- spatial_pathway_test_object()
  expect_warning(
    plot <- PlotSinglePathwaySpatially(
      "unmeasured", object, "slice1", title = "Missing pathway"
    ),
    "No pathway features.*RNA.*scale.data"
  )
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Missing pathway")
  expect_match(ggplot2::ggplot_build(plot)$data[[1]]$label, "No pathway features")
  expect_false("pathway_x" %in% colnames(object[[]]))
})

test_that("spatial pathways with undefined z-scores return a labelled panel", {
  object <- spatial_pathway_test_object()
  expect_warning(
    plot <- PlotSinglePathwaySpatially(
      "constant", object, "slice1", title = "Constant pathway"
    ),
    "No finite pathway z-scores"
  )
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Constant pathway")
  expect_match(ggplot2::ggplot_build(plot)$data[[1]]$label, "No finite pathway z-scores")
})

test_that("valid FOV pathway plots preserve z-scores and render", {
  object <- spatial_pathway_test_object()
  plot <- PlotSinglePathwaySpatially(
    c("varying", "unmeasured"), object, "slice1", title = "Measured pathway"
  )
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Measured pathway")
  expected <- as.numeric(scale(c(1, 3, 2, 5)))
  expect_equal(sort(as.numeric(plot$data$pathway_x)), sort(expected))
  expect_equal(plot$scales$get_scales("fill")$limits, c(-3, 3))
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})

spatial_pathway_test_database <- function() {
  list(
    pathway = data.frame(pathwayRampId = c("p1", "p2"),
                         pathwayName = c("Measured", "Unmeasured")),
    analytehaspathway = data.frame(pathwayRampId = c("p1", "p2"),
                                   rampId = c("varying", "absent"))
  )
}

test_that("spatial pathway batches preserve missing panels and pathway order", {
  object <- spatial_pathway_test_object()
  expect_warning(
    plots <- PlotPathwaysSpatially(
      c("Unmeasured", "Measured"), object, "slice1",
      database = spatial_pathway_test_database()
    ),
    "No pathway features"
  )
  expect_named(plots, c("Unmeasured", "Measured"))
  expect_identical(plots[[1]]$labels$title, "Unmeasured")
  expect_identical(plots[[2]]$labels$title, "Measured")
  expect_s3_class(ggplot2::ggplotGrob(plots[[1]]), "gtable")
})

test_that("unknown spatial pathway names fail before score calculation", {
  expect_error(
    PlotPathwaysSpatially(
      c("Measured", "Typo"), spatial_pathway_test_object(), "slice1",
      database = spatial_pathway_test_database()
    ),
    "not found.*Typo"
  )
})
