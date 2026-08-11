make_alignment_object <- function(x, y, cells = paste0("cell", seq_along(x))) {
  counts <- Matrix::Matrix(matrix(
    seq_len(2L * length(cells)),
    nrow = 2L,
    dimnames = list(c("feature1", "feature2"), cells)
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  centroids <- SeuratObject::CreateCentroids(
    data.frame(x = x, y = y, cell = cells)
  )
  object[["fov"]] <- SeuratObject::CreateFOV(
    coords = list(centroids = centroids),
    type = "centroids",
    assay = SeuratObject::DefaultAssay(object),
    key = "align_"
  )
  object
}


test_that("ApplySpatialAlignment applies SMINT coordinate columns by cell ID", {
  sm <- make_alignment_object(c(0, 1, 2), c(5, 6, 7), c("a", "b", "c"))
  external <- data.frame(
    cell = c("c", "a", "b"),
    x_transformed = c(30, 10, 20),
    y_transformed = c(60, 40, 50)
  )

  result <- ApplySpatialAlignment(
    sm,
    alignment = external,
    return = "result",
    verbose = FALSE
  )

  coordinates <- SeuratObject::GetTissueCoordinates(result$object, image = "fov")
  expect_equal(coordinates$x, c(10, 20, 30))
  expect_equal(coordinates$y, c(40, 50, 60))
  expect_equal(result$coordinates$cell, c("a", "b", "c"))
  expect_equal(result$provenance$backend, "precomputed")
  expect_true("SMINT" %in% names(result$object@tools$spatial_alignment))
})


test_that("ApplySpatialAlignment accepts final coordinate conventions", {
  sm <- make_alignment_object(c(0, 1, 2), c(5, 6, 7))
  external <- data.frame(x_final = c(3, 4, 5), y_final = c(8, 9, 10))

  aligned <- ApplySpatialAlignment(sm, alignment = external, verbose = FALSE)
  coordinates <- SeuratObject::GetTissueCoordinates(aligned, image = "fov")

  expect_equal(coordinates$x, external$x_final)
  expect_equal(coordinates$y, external$y_final)
})


test_that("ApplySpatialAlignment applies homogeneous matrices", {
  sm <- make_alignment_object(c(0, 1, 2), c(5, 6, 7))
  transformation <- matrix(
    c(1, 0, 10, 0, 1, -2, 0, 0, 1),
    nrow = 3L,
    byrow = TRUE
  )

  aligned <- ApplySpatialAlignment(
    sm,
    alignment = list(transformation = list(matrix = transformation)),
    verbose = FALSE
  )
  coordinates <- SeuratObject::GetTissueCoordinates(aligned, image = "fov")

  expect_equal(coordinates$x, c(10, 11, 12))
  expect_equal(coordinates$y, c(3, 4, 5))
})


test_that("ApplySpatialAlignment estimates affine transforms from landmarks", {
  source_x <- c(0, 2, 0, 2)
  source_y <- c(0, 0, 2, 2)
  sm <- make_alignment_object(source_x, source_y)
  st <- make_alignment_object(source_x + 5, source_y - 3)
  source_landmarks <- cbind(x = source_x[1:3], y = source_y[1:3])
  target_landmarks <- source_landmarks + matrix(c(5, -3), nrow = 3, ncol = 2, byrow = TRUE)

  result <- ApplySpatialAlignment(
    sm,
    ST.data = st,
    method = "affine",
    SM.landmarks = source_landmarks,
    ST.landmarks = target_landmarks,
    return = "result",
    verbose = FALSE
  )

  expect_equal(result$coordinates$x_aligned, source_x + 5, tolerance = 1e-10)
  expect_equal(result$coordinates$y_aligned, source_y - 3, tolerance = 1e-10)
  expect_equal(result$diagnostics$landmark_rmse, 0, tolerance = 1e-10)
  expect_lt(
    result$diagnostics$nearest_target_after$median,
    result$diagnostics$nearest_target_before$median
  )
})


test_that("alignment preprocessing reproduces scale and left rotation", {
  coordinates <- data.frame(x = c(1, 2), y = c(3, 4))
  transformed <- .sa_preprocess_coordinates(
    coordinates,
    scale = 10,
    rotate = 90,
    origin = c(0, 0)
  )

  expect_equal(transformed$coordinates$x, c(10, 0), tolerance = 1e-10)
  expect_equal(transformed$coordinates$y, c(10, 20), tolerance = 1e-10)
})


test_that("point-annotator y-x landmarks can be marked as preprocessed", {
  preprocessing <- .sa_preprocess_coordinates(
    data.frame(x = c(1, 2, 1), y = c(3, 3, 4)),
    scale = 10,
    rotate = 90,
    origin = c(0, 0)
  )
  processed_source <- as.matrix(preprocessing$coordinates)
  target <- processed_source + 5

  landmarks <- .sa_prepare_landmarks(
    SM.landmarks = processed_source[, c(2, 1)],
    ST.landmarks = target[, c(2, 1)],
    preprocessing = preprocessing$matrix,
    landmark_order = "yx",
    SM_landmark_space = "preprocessed"
  )

  expect_equal(landmarks$source, processed_source)
  expect_equal(landmarks$target, target)
})


test_that("ApplySpatialAlignment validates incomplete external results", {
  sm <- make_alignment_object(c(0, 1, 2), c(5, 6, 7), c("a", "b", "c"))

  expect_error(
    ApplySpatialAlignment(
      sm,
      alignment = data.frame(cell = c("a", "b"), x = 1:2, y = 3:4),
      verbose = FALSE
    ),
    "missing 1 SM cell identifier",
    fixed = TRUE
  )
  expect_error(
    ApplySpatialAlignment(sm, alignment = data.frame(foo = 1:3), verbose = FALSE),
    "Could not identify aligned x-y columns",
    fixed = TRUE
  )
})


test_that("STalign integration is opt-in for local backend testing", {
  skip_if_not(identical(Sys.getenv("SPAMTP_TEST_STALIGN"), "true"))
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("STalign"))

  source <- expand.grid(x = seq(0, 40, 10), y = seq(0, 40, 10))
  target <- transform(source, x = x + 5, y = y - 3)
  sm <- make_alignment_object(source$x, source$y)
  st <- make_alignment_object(target$x, target$y)

  result <- ApplySpatialAlignment(
    sm,
    ST.data = st,
    method = "lddmm",
    SM.landmarks = source,
    ST.landmarks = target,
    niter = 1,
    lddmm.args = list(a = 10, epV = 20),
    device = "cpu",
    return = "result",
    verbose = FALSE
  )

  expect_equal(result$coordinates$x_aligned, target$x, tolerance = 1e-2)
  expect_equal(result$coordinates$y_aligned, target$y, tolerance = 1e-2)
  expect_null(result$transformation)
  expect_equal(result$backend$engine, "STalign")
  expect_equal(dim(result$backend$affine_yx), c(3, 3))
})
