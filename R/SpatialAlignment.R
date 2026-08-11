#' Apply automatic or precomputed spatial alignment
#'
#' Aligns a moving Spatial Metabolomics (SM) Seurat object to a fixed Spatial
#' Transcriptomics (ST) object. The function can apply coordinates exported by
#' the SMINT workflow, apply a homogeneous transformation matrix, estimate an
#' affine transform from paired landmarks, or run the STalign LDDMM backend used
#' by SMINT through `reticulate`.
#'
#' Python is only required when `method = "lddmm"` and `alignment = NULL`.
#' Precomputed SMINT coordinates and affine landmark alignment are handled
#' entirely in R. The backend summary records STalign's affine component as
#' `affine_yx`; it is not the complete nonlinear LDDMM transformation.
#'
#' @param SM.data A Seurat object containing the moving SM coordinates.
#' @param ST.data An optional Seurat object containing the fixed ST coordinates.
#'   It is required when computing an affine or LDDMM alignment, but is optional
#'   when applying final coordinates supplied through `alignment`.
#' @param alignment Optional precomputed alignment. Accepted inputs are a data
#'   frame or matrix of final coordinates; a 2-by-3 or 3-by-3 homogeneous
#'   transformation matrix; a list containing `aligned_coordinates`,
#'   `coordinates`, `transformation`, or `matrix`; a CSV/TSV/RDS/JSON file; or a
#'   SMINT output directory containing `aligned_coordinates.csv` and optionally
#'   `transformation.json`. Coordinate columns produced by the existing SpaMTP
#'   and SMINT notebooks (`x_transformed`, `x_new`, and `x_final`, with their y
#'   counterparts) are detected automatically.
#' @param method Alignment method used when `alignment` is `NULL`. `"affine"`
#'   estimates a six-degree-of-freedom transformation from paired landmarks in
#'   R. `"lddmm"` runs STalign's landmark-guided affine plus diffeomorphic
#'   registration through Python.
#' @param SM.fov Name of the FOV containing the moving SM centroids. By default,
#'   the first image in `SM.data` is used.
#' @param ST.image Name of the image/FOV containing fixed ST coordinates. By
#'   default, the first image in `ST.data` is used.
#' @param SM.boundary Name of the centroid boundary inside `SM.fov`.
#' @param SM.landmarks,ST.landmarks Matched landmark coordinates as two-column
#'   matrices or data frames in x-y order. At least three non-collinear pairs are
#'   required for `method = "affine"`. Landmarks are optional but strongly
#'   recommended for `method = "lddmm"` when the sections are not already close.
#'   ST landmarks must use the target coordinate system after any
#'   `ST.scale.factor` has been applied.
#' @param landmark.order Coordinate-column order in both landmark inputs.
#'   `"xy"` is the R/SpaMTP convention. Use `"yx"` for row-column arrays saved
#'   by the notebook's point annotator.
#' @param SM.landmark.space Whether SM landmarks use the `"original"` FOV
#'   coordinates or the `"preprocessed"` coordinates after source scaling,
#'   rotation, reflection, translation, and origin adjustment. Point-annotator
#'   landmarks made from the rasterized notebook output are preprocessed.
#' @param coordinate.columns Optional names or integer positions of the x and y
#'   columns in a precomputed coordinate table.
#' @param cell.column Optional name or integer position of the cell/spot ID
#'   column in a precomputed coordinate table. When no IDs are available, rows
#'   must already be in the same order as the moving SM centroids.
#' @param ST.scale.factor Optional fixed-coordinate scale factor. Supply a
#'   numeric value or a scale-factor name such as `"hires"` or `"lowres"` for
#'   Visium images. `NULL` uses full-resolution coordinates.
#' @param source.scale Numeric scalar or x-y vector applied to moving coordinates
#'   before fitting. This reproduces the initial SM scaling step in the SMINT
#'   notebooks.
#' @param source.rotate Counter-clockwise rotation in degrees, around the origin,
#'   applied after scaling.
#' @param source.flip Logical x-y vector indicating reflection around the origin.
#' @param source.translate Numeric x-y translation applied after scaling,
#'   reflection, and rotation.
#' @param source.origin Optional numeric x-y lower bound. Coordinates are shifted
#'   only where necessary so that their minima are at least this value. Use
#'   `c(0, 0)` to reproduce the non-negative coordinate step in the notebook.
#' @param dx Raster spacing passed to `STalign.rasterize()`.
#' @param niter Number of LDDMM optimisation iterations.
#' @param epV Velocity-field gradient step passed to `STalign.LDDMM()`.
#' @param device One of `"auto"`, `"cpu"`, or `"cuda"`. `"auto"` uses CUDA
#'   when Torch reports that it is available.
#' @param seed Random seed set for NumPy and Torch before LDDMM.
#' @param lddmm.args Named list of additional arguments passed to
#'   `STalign.LDDMM()`. Values here override `niter` and `epV`.
#' @param python Optional path to a Python executable containing `STalign` and
#'   `torch`. It must be selected before reticulate initialises Python.
#' @param alignment.name Name used to store alignment provenance in
#'   `SM.data@tools$spatial_alignment`.
#' @param store Logical; store lightweight provenance and diagnostics in the
#'   returned object.
#' @param return Either `"object"` (default) or `"result"`. The latter returns a
#'   list containing the object, coordinate tables, diagnostics, parameters, and
#'   lightweight backend information.
#' @param diagnostics.max.points Maximum number of moving and fixed points used
#'   for nearest-neighbour diagnostics.
#' @param verbose Logical; print progress messages.
#'
#' @return A Seurat object when `return = "object"`; otherwise an object of class
#'   `spamtp_spatial_alignment` containing the aligned object and alignment
#'   details.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Apply coordinates already exported by the SMINT/STalign notebook.
#' aligned_sm <- ApplySpatialAlignment(
#'   SM.data = sm,
#'   ST.data = st,
#'   alignment = "Ven5_z2_transformed_metabolites_coordinates.csv",
#'   coordinate.columns = c("x_transformed", "y_transformed")
#' )
#'
#' # Run the notebook's landmark-guided LDDMM workflow from R.
#' fit <- ApplySpatialAlignment(
#'   SM.data = sm,
#'   ST.data = st,
#'   method = "lddmm",
#'   SM.landmarks = sm_landmarks,
#'   ST.landmarks = st_landmarks,
#'   landmark.order = "yx",             # point-annotator row/column order
#'   SM.landmark.space = "preprocessed", # landmarks selected after rasterization
#'   source.scale = 10,
#'   source.rotate = 90,
#'   source.origin = c(0, 0),
#'   python = "/path/to/STalign_env/bin/python",
#'   return = "result"
#' )
#' fit$diagnostics
#' }
ApplySpatialAlignment <- function(
    SM.data,
    ST.data = NULL,
    alignment = NULL,
    method = c("lddmm", "affine"),
    SM.fov = NULL,
    ST.image = NULL,
    SM.boundary = "centroids",
    SM.landmarks = NULL,
    ST.landmarks = NULL,
    landmark.order = c("xy", "yx"),
    SM.landmark.space = c("original", "preprocessed"),
    coordinate.columns = NULL,
    cell.column = NULL,
    ST.scale.factor = NULL,
    source.scale = 1,
    source.rotate = 0,
    source.flip = c(FALSE, FALSE),
    source.translate = c(0, 0),
    source.origin = NULL,
    dx = 30,
    niter = 1000,
    epV = 200,
    device = c("auto", "cpu", "cuda"),
    seed = 1,
    lddmm.args = list(),
    python = NULL,
    alignment.name = "SMINT",
    store = TRUE,
    return = c("object", "result"),
    diagnostics.max.points = 10000,
    verbose = TRUE
) {
  method <- match.arg(method)
  device <- match.arg(device)
  return <- match.arg(return)
  landmark.order <- match.arg(landmark.order)
  SM.landmark.space <- match.arg(SM.landmark.space)
  if (!is.character(alignment.name) || length(alignment.name) != 1L ||
      is.na(alignment.name) || !nzchar(alignment.name)) {
    stop("`alignment.name` must be one non-empty character string.", call. = FALSE)
  }
  for (argument in c("store", "verbose")) {
    value <- get(argument)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  .sa_validate_seurat(SM.data, "SM.data")
  SM.fov <- .sa_resolve_image(SM.data, SM.fov, "SM.fov")
  source <- .sa_get_coordinates(SM.data, SM.fov)

  preprocessing <- .sa_preprocess_coordinates(
    source[, c("x", "y"), drop = FALSE],
    scale = source.scale,
    rotate = source.rotate,
    flip = source.flip,
    translate = source.translate,
    origin = source.origin
  )
  moving <- preprocessing$coordinates
  moving$cell <- source$cell

  target <- NULL
  if (!is.null(ST.data)) {
    .sa_validate_seurat(ST.data, "ST.data")
    ST.image <- .sa_resolve_image(ST.data, ST.image, "ST.image")
    target <- .sa_get_coordinates(ST.data, ST.image)
    target_scale <- .sa_scale_factor(ST.data, ST.image, ST.scale.factor)
    target[, c("x", "y")] <- target[, c("x", "y"), drop = FALSE] * target_scale
  }

  landmarks <- .sa_prepare_landmarks(
    SM.landmarks = SM.landmarks,
    ST.landmarks = ST.landmarks,
    preprocessing = preprocessing$matrix,
    landmark_order = landmark.order,
    SM_landmark_space = SM.landmark.space
  )

  backend <- "precomputed"
  backend_info <- list()
  transformation <- NULL
  transformed_landmarks <- NULL

  if (!is.null(alignment)) {
    payload <- .sa_read_alignment(alignment)
    coordinates_payload <- .sa_find_coordinates(payload)
    transformation <- .sa_find_transformation(payload)

    if (!is.null(coordinates_payload)) {
      aligned <- .sa_match_coordinates(
        coordinates_payload,
        source_cells = source$cell,
        coordinate_columns = coordinate.columns,
        cell_column = cell.column
      )
    } else if (!is.null(transformation)) {
      aligned_xy <- .sa_apply_transform(moving[, c("x", "y"), drop = FALSE], transformation)
      aligned <- data.frame(cell = source$cell, x = aligned_xy[, 1], y = aligned_xy[, 2])
      if (!is.null(landmarks)) {
        transformed_landmarks <- .sa_apply_transform(landmarks$source, transformation)
      }
    } else {
      stop(
        "`alignment` did not contain aligned coordinates or a supported transformation matrix.",
        call. = FALSE
      )
    }
  } else {
    if (is.null(target)) {
      stop("`ST.data` is required when `alignment` is NULL.", call. = FALSE)
    }

    if (identical(method, "affine")) {
      if (is.null(landmarks)) {
        stop("Paired `SM.landmarks` and `ST.landmarks` are required for affine alignment.", call. = FALSE)
      }
      transformation <- .sa_fit_affine(landmarks$source, landmarks$target)
      aligned_xy <- .sa_apply_transform(moving[, c("x", "y"), drop = FALSE], transformation)
      transformed_landmarks <- .sa_apply_transform(landmarks$source, transformation)
      aligned <- data.frame(cell = source$cell, x = aligned_xy[, 1], y = aligned_xy[, 2])
      backend <- "R-affine"
    } else {
      if (verbose) {
        message("Running SMINT-compatible STalign LDDMM registration ...")
      }
      fit <- .sa_run_stalign(
        moving = moving[, c("x", "y"), drop = FALSE],
        target = target[, c("x", "y"), drop = FALSE],
        landmarks = landmarks,
        dx = dx,
        niter = niter,
        epV = epV,
        device = device,
        seed = seed,
        lddmm.args = lddmm.args,
        python = python
      )
      aligned <- data.frame(
        cell = source$cell,
        x = fit$coordinates[, 1],
        y = fit$coordinates[, 2]
      )
      transformation <- NULL
      transformed_landmarks <- fit$landmarks
      backend_info <- fit$backend
      backend <- "STalign-LDDMM"
    }
  }

  if (any(!is.finite(as.matrix(aligned[, c("x", "y"), drop = FALSE])))) {
    stop("Aligned coordinates contain non-finite values.", call. = FALSE)
  }

  diagnostics <- .sa_alignment_diagnostics(
    original = source[, c("x", "y"), drop = FALSE],
    aligned = aligned[, c("x", "y"), drop = FALSE],
    target = if (is.null(target)) NULL else target[, c("x", "y"), drop = FALSE],
    transformed_landmarks = transformed_landmarks,
    target_landmarks = if (is.null(landmarks)) NULL else landmarks$target,
    max_points = diagnostics.max.points
  )

  output <- .sa_set_fov_coordinates(
    SM.data,
    fov = SM.fov,
    boundary = SM.boundary,
    aligned = aligned
  )

  parameters <- list(
    method = if (is.null(alignment)) method else "precomputed",
    backend = backend,
    SM.fov = SM.fov,
    ST.image = ST.image,
    ST.scale.factor = ST.scale.factor,
    source.scale = source.scale,
    source.rotate = source.rotate,
    source.flip = source.flip,
    source.translate = source.translate,
    source.origin = source.origin,
    landmark.order = landmark.order,
    SM.landmark.space = SM.landmark.space,
    dx = if (identical(backend, "STalign-LDDMM")) dx else NULL,
    niter = if (identical(backend, "STalign-LDDMM")) niter else NULL,
    epV = if (identical(backend, "STalign-LDDMM")) epV else NULL,
    device = if (identical(backend, "STalign-LDDMM")) device else NULL,
    seed = if (identical(backend, "STalign-LDDMM")) seed else NULL,
    lddmm.args = if (identical(backend, "STalign-LDDMM")) lddmm.args else NULL
  )

  provenance <- list(
    name = alignment.name,
    engine = "SpaMTP::ApplySpatialAlignment",
    backend = backend,
    parameters = parameters,
    preprocessing_matrix = preprocessing$matrix,
    transformation = transformation,
    diagnostics = diagnostics,
    backend_info = backend_info
  )

  if (isTRUE(store)) {
    if (is.null(output@tools$spatial_alignment)) {
      output@tools$spatial_alignment <- list()
    }
    output@tools$spatial_alignment[[alignment.name]] <- provenance
  }

  coordinates <- data.frame(
    cell = source$cell,
    x_original = source$x,
    y_original = source$y,
    x_aligned = aligned$x,
    y_aligned = aligned$y,
    stringsAsFactors = FALSE
  )

  result <- structure(
    list(
      object = output,
      coordinates = coordinates,
      target_coordinates = target,
      diagnostics = diagnostics,
      parameters = parameters,
      preprocessing_matrix = preprocessing$matrix,
      transformation = transformation,
      backend = backend_info,
      provenance = provenance
    ),
    class = "spamtp_spatial_alignment"
  )

  if (verbose) {
    message(
      "Aligned ", nrow(aligned), " SM coordinates using ", backend,
      if (!is.null(diagnostics$nearest_target_after$median)) {
        paste0(
          "; median nearest-target distance ",
          signif(diagnostics$nearest_target_before$median, 4), " -> ",
          signif(diagnostics$nearest_target_after$median, 4)
        )
      } else {
        ""
      }
    )
  }

  if (identical(return, "object")) output else result
}


.sa_validate_seurat <- function(object, argument) {
  if (!inherits(object, "Seurat")) {
    stop("`", argument, "` must be a Seurat object.", call. = FALSE)
  }
  invisible(TRUE)
}


.sa_resolve_image <- function(object, image, argument) {
  images <- names(object@images)
  if (!length(images)) {
    stop("`", argument, "` cannot be resolved because the object has no spatial images/FOVs.", call. = FALSE)
  }
  if (is.null(image)) image <- images[[1L]]
  if (length(image) != 1L || is.na(image) || !image %in% images) {
    stop(
      "Invalid `", argument, "`. Available values: ",
      paste(images, collapse = ", "), ".",
      call. = FALSE
    )
  }
  image
}


.sa_get_coordinates <- function(object, image) {
  coordinates <- as.data.frame(SeuratObject::GetTissueCoordinates(object, image = image))
  if (!all(c("x", "y") %in% names(coordinates))) {
    stop("Spatial coordinates for `", image, "` do not contain x and y columns.", call. = FALSE)
  }
  if (!"cell" %in% names(coordinates)) {
    cells <- rownames(coordinates)
    if (is.null(cells) || any(!nzchar(cells))) {
      stop("Spatial coordinates for `", image, "` do not contain cell identifiers.", call. = FALSE)
    }
    coordinates$cell <- cells
  }
  coordinates <- coordinates[, c("x", "y", "cell"), drop = FALSE]
  coordinates$x <- suppressWarnings(as.numeric(coordinates$x))
  coordinates$y <- suppressWarnings(as.numeric(coordinates$y))
  coordinates$cell <- as.character(coordinates$cell)
  if (any(!is.finite(as.matrix(coordinates[, c("x", "y"), drop = FALSE])))) {
    stop("Spatial coordinates for `", image, "` contain non-finite values.", call. = FALSE)
  }
  if (anyDuplicated(coordinates$cell)) {
    stop("Spatial coordinates for `", image, "` contain duplicated cell identifiers.", call. = FALSE)
  }
  rownames(coordinates) <- NULL
  coordinates
}


.sa_scale_factor <- function(object, image, scale_factor) {
  if (is.null(scale_factor)) return(1)
  if (is.numeric(scale_factor)) {
    if (length(scale_factor) != 1L || !is.finite(scale_factor) || scale_factor <= 0) {
      stop("Numeric `ST.scale.factor` must be one finite value greater than zero.", call. = FALSE)
    }
    return(as.numeric(scale_factor))
  }
  if (!is.character(scale_factor) || length(scale_factor) != 1L || is.na(scale_factor)) {
    stop("`ST.scale.factor` must be NULL, numeric, or one scale-factor name.", call. = FALSE)
  }
  spatial_image <- object@images[[image]]
  if (!"scale.factors" %in% methods::slotNames(spatial_image)) {
    stop("Image `", image, "` does not contain named scale factors.", call. = FALSE)
  }
  factors <- spatial_image@scale.factors
  if (!scale_factor %in% names(factors)) {
    stop(
      "Scale factor `", scale_factor, "` was not found. Available values: ",
      paste(names(factors), collapse = ", "), ".",
      call. = FALSE
    )
  }
  value <- as.numeric(factors[[scale_factor]])
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop("Selected ST scale factor is not one finite value greater than zero.", call. = FALSE)
  }
  value
}


.sa_pair <- function(value, argument, mode = c("numeric", "logical")) {
  mode <- match.arg(mode)
  if (length(value) == 1L) value <- rep(value, 2L)
  if (length(value) != 2L || anyNA(value)) {
    stop("`", argument, "` must contain one or two non-missing values.", call. = FALSE)
  }
  if (identical(mode, "numeric")) {
    value <- suppressWarnings(as.numeric(value))
    if (any(!is.finite(value))) {
      stop("`", argument, "` must contain finite numeric values.", call. = FALSE)
    }
  } else {
    if (!is.logical(value)) {
      stop("`", argument, "` must be logical.", call. = FALSE)
    }
  }
  value
}


.sa_apply_transform <- function(points, transformation) {
  points <- as.matrix(points)
  storage.mode(points) <- "double"
  if (ncol(points) != 2L || any(!is.finite(points))) {
    stop("Points must be a finite numeric matrix with two columns.", call. = FALSE)
  }
  transformation <- .sa_matrix(transformation)
  if (identical(dim(transformation), c(2L, 3L))) {
    transformation <- rbind(transformation, c(0, 0, 1))
  }
  if (!identical(dim(transformation), c(3L, 3L))) {
    stop("Transformation must be a 2-by-3 or 3-by-3 numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(transformation))) {
    stop("Transformation must contain only finite numeric values.", call. = FALSE)
  }
  homogeneous <- cbind(points, 1)
  transformed <- homogeneous %*% t(transformation)
  denominator <- transformed[, 3L]
  if (any(!is.finite(denominator)) || any(abs(denominator) < .Machine$double.eps)) {
    stop("Transformation produced invalid homogeneous coordinates.", call. = FALSE)
  }
  transformed[, 1:2, drop = FALSE] / denominator
}


.sa_preprocess_coordinates <- function(
    coordinates,
    scale = 1,
    rotate = 0,
    flip = c(FALSE, FALSE),
    translate = c(0, 0),
    origin = NULL
) {
  scale <- .sa_pair(scale, "source.scale")
  if (any(scale == 0)) stop("`source.scale` cannot contain zero.", call. = FALSE)
  flip <- .sa_pair(flip, "source.flip", mode = "logical")
  translate <- .sa_pair(translate, "source.translate")
  if (!is.numeric(rotate) || length(rotate) != 1L || !is.finite(rotate)) {
    stop("`source.rotate` must be one finite number.", call. = FALSE)
  }

  angle <- rotate * pi / 180
  scaling <- diag(c(scale * ifelse(flip, -1, 1), 1))
  rotation <- matrix(
    c(cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0, 0, 0, 1),
    nrow = 3L
  )
  translation <- diag(3L)
  translation[1:2, 3L] <- translate
  transformation <- translation %*% rotation %*% scaling
  transformed <- .sa_apply_transform(coordinates, transformation)

  if (!is.null(origin)) {
    origin <- .sa_pair(origin, "source.origin")
    shift <- pmax(origin - apply(transformed, 2L, min), 0)
    origin_transform <- diag(3L)
    origin_transform[1:2, 3L] <- shift
    transformation <- origin_transform %*% transformation
    transformed <- .sa_apply_transform(coordinates, transformation)
  }

  colnames(transformed) <- c("x", "y")
  list(coordinates = as.data.frame(transformed), matrix = transformation)
}


.sa_prepare_landmarks <- function(
    SM.landmarks,
    ST.landmarks,
    preprocessing,
    landmark_order,
    SM_landmark_space
) {
  if (is.null(SM.landmarks) && is.null(ST.landmarks)) return(NULL)
  if (is.null(SM.landmarks) || is.null(ST.landmarks)) {
    stop("`SM.landmarks` and `ST.landmarks` must be supplied together.", call. = FALSE)
  }
  source <- .sa_landmark_matrix(SM.landmarks, "SM.landmarks")
  target <- .sa_landmark_matrix(ST.landmarks, "ST.landmarks")
  if (nrow(source) != nrow(target)) {
    stop("`SM.landmarks` and `ST.landmarks` must contain the same number of rows.", call. = FALSE)
  }
  if (identical(landmark_order, "yx")) {
    source <- source[, c(2L, 1L), drop = FALSE]
    target <- target[, c(2L, 1L), drop = FALSE]
  }
  if (identical(SM_landmark_space, "original")) {
    source <- .sa_apply_transform(source, preprocessing)
  }
  colnames(source) <- colnames(target) <- c("x", "y")
  list(source = source, target = target)
}


.sa_landmark_matrix <- function(x, argument) {
  if (inherits(x, "python.builtin.object")) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("Package `reticulate` is required to convert Python landmarks.", call. = FALSE)
    }
    x <- reticulate::py_to_r(x)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (ncol(x) != 2L || nrow(x) < 1L || any(!is.finite(x))) {
    stop("`", argument, "` must be a finite numeric matrix with two columns.", call. = FALSE)
  }
  colnames(x) <- c("x", "y")
  x
}


.sa_fit_affine <- function(source, target) {
  if (nrow(source) < 3L) {
    stop("Affine alignment requires at least three paired landmarks.", call. = FALSE)
  }
  design <- cbind(source, 1)
  if (qr(design)$rank < 3L) {
    stop("Affine landmarks must contain at least three non-collinear source points.", call. = FALSE)
  }
  coefficients <- qr.solve(design, target)
  rbind(t(coefficients), c(0, 0, 1))
}


.sa_read_alignment <- function(alignment) {
  if (inherits(alignment, "python.builtin.object")) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("Package `reticulate` is required to convert a Python alignment object.", call. = FALSE)
    }
    alignment <- reticulate::py_to_r(alignment)
  }
  if (!is.character(alignment) || length(alignment) != 1L) return(alignment)
  if (!file.exists(alignment)) {
    stop("Alignment path does not exist: ", alignment, call. = FALSE)
  }
  if (dir.exists(alignment)) {
    coordinate_file <- file.path(alignment, "aligned_coordinates.csv")
    transform_file <- file.path(alignment, "transformation.json")
    if (!file.exists(coordinate_file) && !file.exists(transform_file)) {
      stop(
        "SMINT alignment directory must contain aligned_coordinates.csv or transformation.json.",
        call. = FALSE
      )
    }
    return(list(
      aligned_coordinates = if (file.exists(coordinate_file)) {
        utils::read.csv(coordinate_file, check.names = FALSE)
      } else {
        NULL
      },
      transformation = if (file.exists(transform_file)) {
        jsonlite::fromJSON(transform_file, simplifyVector = TRUE)
      } else {
        NULL
      }
    ))
  }

  extension <- tolower(sub("^.*[.]", "", basename(alignment)))
  switch(
    extension,
    csv = utils::read.csv(alignment, check.names = FALSE),
    tsv = utils::read.delim(alignment, check.names = FALSE),
    txt = utils::read.delim(alignment, check.names = FALSE),
    rds = readRDS(alignment),
    json = jsonlite::fromJSON(alignment, simplifyVector = TRUE),
    stop("Unsupported alignment file extension: .", extension, call. = FALSE)
  )
}


.sa_find_coordinates <- function(payload) {
  if (is.data.frame(payload) || (is.matrix(payload) && !.sa_is_transform(payload))) {
    return(payload)
  }
  if (!is.list(payload)) return(NULL)
  for (name in c("aligned_coordinates", "coordinates", "aligned.coords", "aligned")) {
    if (!is.null(payload[[name]])) return(payload[[name]])
  }
  NULL
}


.sa_find_transformation <- function(payload) {
  if (.sa_is_transform(payload)) return(.sa_matrix(payload))
  if (!is.list(payload)) return(NULL)
  candidates <- list(payload$transformation, payload$transform, payload$matrix)
  for (candidate in candidates) {
    if (is.null(candidate)) next
    if (is.list(candidate) && !is.null(candidate$matrix)) candidate <- candidate$matrix
    if (.sa_is_transform(candidate)) return(.sa_matrix(candidate))
  }
  NULL
}


.sa_matrix <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(x$matrix)) x <- x$matrix
    if (is.list(x) && length(x) && all(vapply(x, is.atomic, logical(1)))) {
      x <- do.call(rbind, x)
    }
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}


.sa_is_transform <- function(x) {
  if (is.data.frame(x)) return(FALSE)
  candidate <- suppressWarnings(try(.sa_matrix(x), silent = TRUE))
  !inherits(candidate, "try-error") &&
    (identical(dim(candidate), c(2L, 3L)) || identical(dim(candidate), c(3L, 3L))) &&
    all(is.finite(candidate))
}


.sa_match_coordinates <- function(
    coordinates,
    source_cells,
    coordinate_columns = NULL,
    cell_column = NULL
) {
  coordinates <- as.data.frame(coordinates, check.names = FALSE)
  if (!nrow(coordinates)) stop("Precomputed alignment contains no coordinates.", call. = FALSE)

  if (is.null(coordinate_columns)) {
    pairs <- list(
      c("x_final", "y_final"),
      c("x_transformed", "y_transformed"),
      c("x_new", "y_new"),
      c("x_aligned", "y_aligned"),
      c("aligned_x", "aligned_y"),
      c("x", "y")
    )
    coordinate_columns <- NULL
    for (pair in pairs) {
      if (all(pair %in% names(coordinates))) {
        coordinate_columns <- pair
        break
      }
    }
    if (is.null(coordinate_columns) && ncol(coordinates) == 2L) {
      coordinate_columns <- 1:2
    }
  }
  if (length(coordinate_columns) != 2L ||
      (is.character(coordinate_columns) && !all(coordinate_columns %in% names(coordinates)))) {
    stop(
      "Could not identify aligned x-y columns. Supply `coordinate.columns = c(x, y)`.",
      call. = FALSE
    )
  }
  if (is.numeric(coordinate_columns) &&
      any(coordinate_columns < 1 | coordinate_columns > ncol(coordinates))) {
    stop("Numeric `coordinate.columns` are outside the coordinate table.", call. = FALSE)
  }

  xy <- coordinates[, coordinate_columns, drop = FALSE]
  xy[] <- lapply(xy, function(value) suppressWarnings(as.numeric(value)))
  if (any(!is.finite(as.matrix(xy)))) {
    stop("Precomputed aligned coordinates contain non-finite values.", call. = FALSE)
  }

  if (is.null(cell_column)) {
    candidates <- c("cell", "barcode", "spot", "id", "Cell", "Barcode")
    matches <- candidates[candidates %in% names(coordinates)]
    cell_column <- if (length(matches)) matches[[1L]] else NULL
  }

  ids <- NULL
  if (!is.null(cell_column)) {
    if (length(cell_column) != 1L ||
        (is.character(cell_column) && !cell_column %in% names(coordinates))) {
      stop("Invalid `cell.column`.", call. = FALSE)
    }
    if (is.numeric(cell_column) &&
        (cell_column < 1 || cell_column > ncol(coordinates))) {
      stop("Numeric `cell.column` is outside the coordinate table.", call. = FALSE)
    }
    ids <- as.character(coordinates[[cell_column]])
  } else {
    candidate_rownames <- rownames(coordinates)
    default_rownames <- identical(candidate_rownames, as.character(seq_len(nrow(coordinates))))
    if (!is.null(candidate_rownames) && !default_rownames && all(source_cells %in% candidate_rownames)) {
      ids <- candidate_rownames
    }
  }

  if (is.null(ids)) {
    if (nrow(xy) != length(source_cells)) {
      stop(
        "Precomputed coordinates have no cell IDs and their row count does not match the SM FOV.",
        call. = FALSE
      )
    }
    ids <- source_cells
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Precomputed alignment contains missing or duplicated cell identifiers.", call. = FALSE)
  }
  index <- match(source_cells, ids)
  if (anyNA(index)) {
    stop(
      "Precomputed alignment is missing ", sum(is.na(index)),
      " SM cell identifier(s).",
      call. = FALSE
    )
  }
  data.frame(
    cell = source_cells,
    x = xy[[1L]][index],
    y = xy[[2L]][index],
    stringsAsFactors = FALSE
  )
}


.sa_set_fov_coordinates <- function(object, fov, boundary, aligned) {
  spatial_image <- object[[fov]]
  if (!inherits(spatial_image, "FOV")) {
    stop(
      "`ApplySpatialAlignment()` currently updates Seurat FOV centroids; `", fov,
      "` has class ", paste(class(spatial_image), collapse = "/"), ".",
      call. = FALSE
    )
  }
  boundary_names <- names(spatial_image@boundaries)
  if (!boundary %in% boundary_names) {
    stop(
      "Boundary `", boundary, "` was not found in FOV `", fov,
      "`. Available values: ", paste(boundary_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  centroids <- spatial_image[[boundary]]
  if (!inherits(centroids, "Centroids")) {
    stop("Boundary `", boundary, "` is not a Centroids object.", call. = FALSE)
  }
  cells <- as.character(SeuratObject::Cells(centroids))
  index <- match(cells, aligned$cell)
  if (anyNA(index)) {
    stop("Aligned coordinates are missing cells stored in the FOV boundary.", call. = FALSE)
  }
  coordinates <- as.matrix(aligned[index, c("x", "y"), drop = FALSE])
  storage.mode(coordinates) <- "double"
  colnames(coordinates) <- c("x", "y")
  rownames(coordinates) <- NULL
  centroids@coords <- coordinates
  methods::validObject(centroids)
  spatial_image[[boundary]] <- centroids
  methods::validObject(spatial_image)
  object[[fov]] <- spatial_image
  object
}


.sa_run_stalign <- function(
    moving,
    target,
    landmarks,
    dx,
    niter,
    epV,
    device,
    seed,
    lddmm.args,
    python
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package `reticulate` is required for `method = \"lddmm\"`. ",
      "Install it or apply precomputed SMINT coordinates.",
      call. = FALSE
    )
  }
  for (entry in list(dx = dx, niter = niter, epV = epV, seed = seed)) {
    if (!is.numeric(entry) || length(entry) != 1L || !is.finite(entry)) {
      stop("`dx`, `niter`, `epV`, and `seed` must be finite scalar values.", call. = FALSE)
    }
  }
  if (dx <= 0 || niter < 1 || epV <= 0) {
    stop("`dx`, `niter`, and `epV` must be greater than zero.", call. = FALSE)
  }
  if (niter != floor(niter) || seed != floor(seed)) {
    stop("`niter` and `seed` must be whole numbers.", call. = FALSE)
  }
  if (!is.list(lddmm.args) || is.null(names(lddmm.args)) && length(lddmm.args)) {
    stop("`lddmm.args` must be a named list.", call. = FALSE)
  }
  if (length(lddmm.args) && any(!nzchar(names(lddmm.args)))) {
    stop("Every entry in `lddmm.args` must have a name.", call. = FALSE)
  }
  protected <- c("xI", "I", "xJ", "J", "pointsI", "pointsJ", "L", "T", "A", "v", "xv", "device")
  if (any(names(lddmm.args) %in% protected)) {
    stop(
      "`lddmm.args` cannot override internal arguments: ",
      paste(intersect(names(lddmm.args), protected), collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.null(python)) {
    if (!is.character(python) || length(python) != 1L || is.na(python) ||
        !file.exists(python)) {
      stop("`python` must point to an existing Python executable.", call. = FALSE)
    }
    reticulate::use_python(python, required = TRUE)
  }
  if (!reticulate::py_module_available("STalign")) {
    stop(
      "Python module `STalign` is unavailable. Install STalign in the Python ",
      "environment selected by reticulate, or supply precomputed coordinates.",
      call. = FALSE
    )
  }
  if (!reticulate::py_module_available("torch")) {
    stop("Python module `torch` is required by STalign.", call. = FALSE)
  }

  stalign_module <- reticulate::import("STalign", convert = FALSE)
  stalign <- stalign_module$STalign
  np <- reticulate::import("numpy", convert = FALSE)
  torch <- reticulate::import("torch", convert = FALSE)

  np$random$seed(as.integer(seed))
  torch$manual_seed(as.integer(seed))
  cuda_available <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  selected_device <- switch(
    device,
    auto = if (cuda_available) "cuda:0" else "cpu",
    cpu = "cpu",
    cuda = "cuda:0"
  )
  if (identical(device, "cuda") && !cuda_available) {
    stop("`device = \"cuda\"` was requested but Torch reports no available CUDA device.", call. = FALSE)
  }
  if (cuda_available) torch$cuda$manual_seed_all(as.integer(seed))

  moving_x <- np$array(moving$x, dtype = "float64")
  moving_y <- np$array(moving$y, dtype = "float64")
  target_x <- np$array(target$x, dtype = "float64")
  target_y <- np$array(target$y, dtype = "float64")

  moving_raster <- stalign$rasterize(moving_x, moving_y, dx = as.numeric(dx))
  target_raster <- stalign$rasterize(target_x, target_y, dx = as.numeric(dx))
  xI <- reticulate::py_get_item(moving_raster, 0L)
  yI <- reticulate::py_get_item(moving_raster, 1L)
  I <- reticulate::py_get_item(moving_raster, 2L)
  xJ <- reticulate::py_get_item(target_raster, 0L)
  yJ <- reticulate::py_get_item(target_raster, 1L)
  J <- reticulate::py_get_item(target_raster, 2L)

  moving_yx <- np$stack(reticulate::tuple(moving_y, moving_x), axis = 1L)
  pointsI <- NULL
  pointsJ <- NULL
  initial_L <- NULL
  initial_T <- NULL
  if (!is.null(landmarks)) {
    pointsI <- np$array(landmarks$source[, c(2L, 1L), drop = FALSE], dtype = "float64")
    pointsJ <- np$array(landmarks$target[, c(2L, 1L), drop = FALSE], dtype = "float64")
    initial <- stalign$L_T_from_points(pointsI, pointsJ)
    initial_L <- reticulate::py_get_item(initial, 0L)
    initial_T <- reticulate::py_get_item(initial, 1L)
  }

  arguments <- list(
    reticulate::tuple(yI, xI),
    I,
    reticulate::tuple(yJ, xJ),
    J,
    pointsI = pointsI,
    pointsJ = pointsJ,
    L = initial_L,
    T = initial_T,
    niter = as.integer(niter),
    epV = as.numeric(epV),
    device = selected_device
  )
  arguments[names(lddmm.args)] <- lddmm.args

  fit <- tryCatch(
    do.call(stalign$LDDMM, arguments),
    error = function(error) {
      stop(
        "STalign LDDMM failed: ", conditionMessage(error),
        ". Check coordinate scale, landmarks, `dx`, and `lddmm.args`.",
        call. = FALSE
      )
    }
  )

  transformed <- stalign$transform_points_source_to_target(
    fit$xv, fit$v, fit$A, moving_yx
  )
  transformed <- reticulate::py_to_r(transformed$detach()$cpu()$numpy())
  aligned <- transformed[, c(2L, 1L), drop = FALSE]
  colnames(aligned) <- c("x", "y")

  transformed_landmarks <- NULL
  if (!is.null(pointsI)) {
    transformed_landmarks <- stalign$transform_points_source_to_target(
      fit$xv, fit$v, fit$A, pointsI
    )
    transformed_landmarks <- reticulate::py_to_r(
      transformed_landmarks$detach()$cpu()$numpy()
    )[, c(2L, 1L), drop = FALSE]
    colnames(transformed_landmarks) <- c("x", "y")
  }

  affine <- reticulate::py_to_r(fit$A$detach()$cpu()$numpy())
  velocity_shape <- unlist(reticulate::py_to_r(fit$v$shape), use.names = FALSE)
  python_version <- reticulate::py_config()$version

  try({
    pyplot <- reticulate::import("matplotlib.pyplot", convert = FALSE)
    pyplot$close("all")
  }, silent = TRUE)

  list(
    coordinates = aligned,
    landmarks = transformed_landmarks,
    backend = list(
      engine = "STalign",
      python = as.character(python_version),
      device = selected_device,
      cuda_available = cuda_available,
      affine_yx = affine,
      velocity_shape = velocity_shape
    )
  )
}


.sa_alignment_diagnostics <- function(
    original,
    aligned,
    target,
    transformed_landmarks,
    target_landmarks,
    max_points
) {
  if (!is.numeric(max_points) || length(max_points) != 1L ||
      !is.finite(max_points) || max_points < 2) {
    stop("`diagnostics.max.points` must be one finite number >= 2.", call. = FALSE)
  }
  bounds <- function(x) {
    list(
      x = unname(range(x[, 1L])),
      y = unname(range(x[, 2L]))
    )
  }
  output <- list(
    original_bounds = bounds(original),
    aligned_bounds = bounds(aligned),
    target_bounds = if (is.null(target)) NULL else bounds(target),
    nearest_target_before = NULL,
    nearest_target_after = NULL,
    landmark_rmse = NULL
  )
  if (!is.null(target)) {
    sample_rows <- function(n) {
      if (n <= max_points) seq_len(n) else unique(round(seq(1, n, length.out = max_points)))
    }
    target_sample <- as.matrix(target[sample_rows(nrow(target)), , drop = FALSE])
    before_sample <- as.matrix(original[sample_rows(nrow(original)), , drop = FALSE])
    after_sample <- as.matrix(aligned[sample_rows(nrow(aligned)), , drop = FALSE])
    summarize_distance <- function(query) {
      distances <- FNN::get.knnx(target_sample, query, k = 1L)$nn.dist[, 1L]
      list(
        median = unname(stats::median(distances)),
        mean = unname(mean(distances)),
        p95 = unname(stats::quantile(distances, 0.95, names = FALSE))
      )
    }
    output$nearest_target_before <- summarize_distance(before_sample)
    output$nearest_target_after <- summarize_distance(after_sample)
  }
  if (!is.null(transformed_landmarks) && !is.null(target_landmarks)) {
    output$landmark_rmse <- sqrt(mean((transformed_landmarks - target_landmarks)^2))
  }
  output
}
