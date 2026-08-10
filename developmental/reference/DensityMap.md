# Generates interactive 3D spatial density plot for m/z values

This function generates an interactive HTML file that can be used to
generate spatial density kernals based on the intensity of any selected
m/z value.

## Usage

``` r
DensityMap(object, assay = "SPM", slot = "counts", folder = getwd(), ...)
```

## Arguments

- object:

  SpaMTP Seurat object contains spatial metabolomics data

- assay:

  Character string defining the Seurat assay that contains intensity
  values to plot (default = "SPM").

- slot:

  Character string defining the assay slot to extract the intensity
  values from (default = "counts").

- folder:

  Character string defining the directory path to save the density map
  HTML file to (default = getwd()).

- ...:

  Arguments passed to SpaMTP::AnnotateSeuratMALDI()

## Value

Return a html contains the annotation/average intensity of
peakbins/density of distribution of each peak

## Examples

``` r
# DensityMap(SpaMTP.obj)
```
