# Perform Dimensionality Reduction using Graph-Regularised PCA on Spatial Data

Computes a graph-regularised PCA using spatial coordinates and scaled
expression data. A k-nearest neighbour (k-NN) graph is computed using
spatial locations and used to regularise the PCA decomposition via a
graph Laplacian. The result are stored in the `@reductions` section of
the returned SpaMTP Seurat object.

## Usage

``` r
RunSpatialGraphPCA(
  data,
  n_components = 50,
  assay = "Spatial",
  slot = "scale.data",
  image = NULL,
  platform = "Visium",
  lambda = 0.5,
  n_neighbors = NULL,
  include_self = FALSE,
  alg = "kd_tree",
  fast = TRUE,
  graph_name = "SpatialKNN",
  reduction_name = "SpatialPCA",
  verbose = TRUE
)
```

## Arguments

- data:

  A SpaMTP Seurat object containing spatial data (feature data and
  spatial coordinates).

- n_components:

  Integer specifying the number of principal components to compute
  (default = 50).

- assay:

  Character string defining the name of the assay to use data from
  (default = "Spatial").

- slot:

  Character string defining the name of the slot to extract scaled data
  from (default = "scale.data").

- image:

  Character string matching the name of the image to use for extracting
  spatial coordinates (default = NULL).

- platform:

  Character string matching either `"Visium"` or `"ST"` to determine how
  the k-NN graph is constructed. If "Visium" k-nns will handle the
  hexagon spot arrangement, including setting `n_neighbors` = 6, else
  "ST" assignment will set `n_neighbors` = 4 unless a value is
  specifically provided (default = "Visium").

- lambda:

  Numeric value defining the regularisation parameter that controls the
  influence of the graph Laplacian (default = 0.5).

- n_neighbors:

  Integer value specifying the number of spatial neighbours to use. If
  `NULL`, will default of 6 for "Visium" data and 4 for "ST" platforms
  (default = NULL).

- include_self:

  Boolean logical value indicating whether to include self-connections
  in the graph (default = FALSE).

- alg:

  Character string specifying the algorithm to use for nearest neighbour
  search (passed to
  [`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html)) (default
  ="kd_tree").

- fast:

  Boolean logical value stating whether to use fast approximate
  eigendecomposition via
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html).
  For large datasets this is recommended (default =TRUE).

- graph_name:

  Character string of the name to use for storing the computed spatial
  graph in `@graphs`(default ="SpatialKNN").

- reduction_name:

  Character string stating the name of the dimensionality reduction
  stored in `@reductions`(default ="SpatialPCA").

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A SpaMTP Seurat object with a new graph stored in `@graphs` and
spatially-aware PCA reduction values stored in `@reductions`.

## Details

Note: This method has been adapted from
[GraphPCA](https://doi.org/10.1186/s13059-024-03429-x) python package.

## Examples

``` r
# spamtp_obj <- RunSpatialGraphPCA(spamtp_obj, platform = "Visium")
```
