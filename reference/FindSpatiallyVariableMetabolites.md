# Find Spatially Variable Metabolites

Finds metabolites that display strong spatial patterns using MoransI.
Each m/z value is ranked by MoransI score and the results are stored in
the SpaMTP Seurat Object feature metadata.

## Usage

``` r
FindSpatiallyVariableMetabolites(
  object,
  assay = "SPM",
  slot = "counts",
  image = "slice1",
  nfeatures = 2000,
  max_spots = 5000,
  seed = 1,
  verbose = TRUE
)
```

## Arguments

- object:

  SpaMTP Seurat class object contating the intensity values for each m/z

- assay:

  Character string indicating which Seurat object assay to pull data
  form (default = "SPM").

- slot:

  Character string indicating the assay slot to use to pull expression
  values form (default = "counts").

- image:

  Character string defining the image to extract the tissue coordinates
  from (defualt = "slice1").

- nfeatures:

  Numeric values defining the top number of features to mark as the top
  spatially variable (default = 2000).

- max_spots:

  Maximum number of spots used to construct Seurat's dense Moran's I
  distance matrix. When the object is larger, spots are sampled
  deterministically using `seed`. Set to `NULL` to use every spot
  (default = 5000).

- seed:

  Integer random seed used only when `max_spots` triggers subsampling
  (default = 1).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

Returns SpaMTP object containing the MoransI pvalue and rank stored in
the assay feature meta.data

## Examples

``` r
# SpaMTP.obj <- FindSpatiallyVariableMetabolites(SpaMTP)
```
