# Generates PCA analysis results for a SpaMTP Seurat Object

This function run PCA analaysis on a SpaMTP Seurat Object. The user can
provide a bin/resolution size to increase the bin size and reduce the
dimensionality/noise of the SM dataset prior to calculating PCAs.

## Usage

``` r
RunMetabolicPCA(
  SpaMTP,
  npcs = 30,
  variance_explained_threshold = 0.9,
  assay = "SPM",
  slot = "counts",
  show_variance_plot = FALSE,
  bin_resolution = NULL,
  resolution_units = "ppm",
  bin_method = "sum",
  reduction.name = "pca",
  verbose = TRUE
)
```

## Arguments

- SpaMTP:

  SpaMTP Seurat class object that contains spatial metabolic
  information.

- npcs:

  is an integer value to indicated preferred number of PCs to retain
  (default = 30).

- variance_explained_threshold:

  Numeric value defining the explained variance threshold (default =
  0.9).

- assay:

  Character string defining the SpaMTP assay to extract intensity values
  from (default = "SPM").

- slot:

  Character string defining the assay slot containing the intensity
  values (default = "counts").

- show_variance_plot:

  Boolean indicating weather to display the variance plot output by this
  analysis (default = FALSE).

- bin_resolution:

  Numeric value defining the resolution to use for binning m/z peaks. If
  set to `NULL`, no binning will be performed (default = NULL).

- resolution_units:

  Character string specifying the units of the `bin_resolution`. Either
  'ppm' or 'mz' can be provided. `bin_resolution` must be provided for
  this parameter to be implemented (default = "ppm").

- bin_method:

  Character string defining the method to use for binning respective m/z
  peaks that fall within a bin. Options for this parameter can be one of
  "sum", "mean", "max" or "min". `bin_resolution` must be provided for
  this parameter to be implemented (default = "sum").

- reduction.name:

  Character string indicating the name associated with the PCA results
  stored in the output SpaMTP Seurat object (default = "pca").

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

SpaMTP object with pca results stored in the

## Examples

``` r
## For running PCA on un-adjusted peak bin sizes
# spamtp_obj <- RunMetabolicPCA(spamtp_obj, npcs = 50)

## For running PCA on increased peak bin sizes
# spamtp_obj <- RunMetabolicPCA(spamtp_obj, npcs = 50, bin_resolution = 50)
```
