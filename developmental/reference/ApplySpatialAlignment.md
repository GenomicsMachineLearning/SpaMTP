# Apply automatic or precomputed spatial alignment

Aligns a moving Spatial Metabolomics (SM) Seurat object to a fixed
Spatial Transcriptomics (ST) object. The function can apply coordinates
exported by the SMINT workflow, apply a homogeneous transformation
matrix, estimate an affine transform from paired landmarks, or run the
STalign LDDMM backend used by SMINT through `reticulate`.

## Usage

``` r
ApplySpatialAlignment(
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
)
```

## Arguments

- SM.data:

  A Seurat object containing the moving SM coordinates.

- ST.data:

  An optional Seurat object containing the fixed ST coordinates. It is
  required when computing an affine or LDDMM alignment, but is optional
  when applying final coordinates supplied through `alignment`.

- alignment:

  Optional precomputed alignment. Accepted inputs are a data frame or
  matrix of final coordinates; a 2-by-3 or 3-by-3 homogeneous
  transformation matrix; a list containing `aligned_coordinates`,
  `coordinates`, `transformation`, or `matrix`; a CSV/TSV/RDS/JSON file;
  or a SMINT output directory containing `aligned_coordinates.csv` and
  optionally `transformation.json`. Coordinate columns produced by the
  existing SpaMTP and SMINT notebooks (`x_transformed`, `x_new`, and
  `x_final`, with their y counterparts) are detected automatically.

- method:

  Alignment method used when `alignment` is `NULL`. `"affine"` estimates
  a six-degree-of-freedom transformation from paired landmarks in R.
  `"lddmm"` runs STalign's landmark-guided affine plus diffeomorphic
  registration through Python.

- SM.fov:

  Name of the FOV containing the moving SM centroids. By default, the
  first image in `SM.data` is used.

- ST.image:

  Name of the image/FOV containing fixed ST coordinates. By default, the
  first image in `ST.data` is used.

- SM.boundary:

  Name of the centroid boundary inside `SM.fov`.

- SM.landmarks, ST.landmarks:

  Matched landmark coordinates as two-column matrices or data frames in
  x-y order. At least three non-collinear pairs are required for
  `method = "affine"`. Landmarks are optional but strongly recommended
  for `method = "lddmm"` when sections are not already close. ST
  landmarks must use the target coordinate system after any
  `ST.scale.factor` has been applied.

- landmark.order:

  Coordinate-column order in both landmark inputs. `"xy"` is the
  R/SpaMTP convention. Use `"yx"` for row-column arrays saved by the
  notebook's point annotator.

- SM.landmark.space:

  Whether SM landmarks use the `"original"` FOV coordinates or the
  `"preprocessed"` coordinates after source scaling, rotation,
  reflection, translation, and origin adjustment. Point-annotator
  landmarks made from the rasterized notebook output are preprocessed.

- coordinate.columns:

  Optional names or integer positions of the x and y columns in a
  precomputed coordinate table.

- cell.column:

  Optional name or integer position of the cell/spot ID column in a
  precomputed coordinate table. When no IDs are available, rows must
  already be in the same order as the moving SM centroids.

- ST.scale.factor:

  Optional fixed-coordinate scale factor. Supply a numeric value or a
  scale-factor name such as `"hires"` or `"lowres"` for Visium images.
  `NULL` uses full-resolution coordinates.

- source.scale:

  Numeric scalar or x-y vector applied to moving coordinates before
  fitting. This reproduces the initial SM scaling step in the SMINT
  notebooks.

- source.rotate:

  Counter-clockwise rotation in degrees, around the origin, applied
  after scaling.

- source.flip:

  Logical x-y vector indicating reflection around the origin.

- source.translate:

  Numeric x-y translation applied after scaling, reflection, and
  rotation.

- source.origin:

  Optional numeric x-y lower bound. Coordinates are shifted only where
  necessary so that their minima are at least this value. Use `c(0, 0)`
  to reproduce the non-negative coordinate step in the notebook.

- dx:

  Raster spacing passed to `STalign.rasterize()`.

- niter:

  Number of LDDMM optimisation iterations.

- epV:

  Velocity-field gradient step passed to `STalign.LDDMM()`.

- device:

  One of `"auto"`, `"cpu"`, or `"cuda"`. `"auto"` uses CUDA when Torch
  reports that it is available.

- seed:

  Random seed set for NumPy and Torch before LDDMM.

- lddmm.args:

  Named list of additional arguments passed to `STalign.LDDMM()`. Values
  here override `niter` and `epV`.

- python:

  Optional path to a Python executable containing `STalign` and `torch`.
  It must be selected before reticulate initialises Python.

- alignment.name:

  Name used to store alignment provenance in
  `SM.data@tools$spatial_alignment`.

- store:

  Logical; store lightweight provenance and diagnostics in the returned
  object.

- return:

  Either `"object"` (default) or `"result"`. The latter returns a list
  containing the object, coordinate tables, diagnostics, parameters, and
  lightweight backend information.

- diagnostics.max.points:

  Maximum number of moving and fixed points used for nearest-neighbour
  diagnostics.

- verbose:

  Logical; print progress messages.

## Value

A Seurat object when `return = "object"`; otherwise an object of class
`spamtp_spatial_alignment` containing the aligned object and alignment
details.

## Details

Python is only required when `method = "lddmm"` and `alignment = NULL`.
Precomputed SMINT coordinates and affine landmark alignment are handled
entirely in R. Lightweight provenance and quality diagnostics can be
stored in the returned object's `@tools$spatial_alignment` entry. The
backend summary records STalign's affine component as `affine_yx`; it is
not the complete nonlinear LDDMM transformation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Apply coordinates already exported by the SMINT/STalign notebook.
aligned_sm <- ApplySpatialAlignment(
  SM.data = sm,
  ST.data = st,
  alignment = "Ven5_z2_transformed_metabolites_coordinates.csv",
  coordinate.columns = c("x_transformed", "y_transformed")
)

# Run the notebook's landmark-guided LDDMM workflow from R.
fit <- ApplySpatialAlignment(
  SM.data = sm,
  ST.data = st,
  method = "lddmm",
  SM.landmarks = sm_landmarks,
  ST.landmarks = st_landmarks,
  landmark.order = "yx",             # point-annotator row/column order
  SM.landmark.space = "preprocessed", # selected after rasterization
  source.scale = 10,
  source.rotate = 90,
  source.origin = c(0, 0),
  python = "/path/to/STalign_env/bin/python",
  return = "result"
)
fit$diagnostics
} # }
```
