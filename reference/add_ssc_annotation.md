# Adds Cardinal ssc segmentation annotation to m/z count data object

Add Spatial Shrunked Centroid (ssc) results from a non-intensity based
Cardinal Object into the pixelData slot of a provided Cardinal Object
containing intensity values for each m/z (e.g. MSImagingExperiment).
This function is additional functionality that can only be run on
`Cardinal` Objects ONLY.

## Usage

``` r
add_ssc_annotation(data, data_ssc, resolution)
```

## Arguments

- data:

  A Cardinal Object containing the raw/binned m/z count data.

- data_ssc:

  A Cardinal Object containing the ssc segmentation results. Note:
  Cardinal's spatialShrunkenCentroids() must be run to generate this
  object.

- resolution:

  An integer (\< v3.6.0) or string (\>= v3.6.0) defining the ssc
  segmentation resolution to add to the raw Cardinal object.

## Value

A Cardinal Object containing the raw m/z counts and the relative ssc
segmentation annotations per pixel.

## Examples

``` r
# ssc_data <- Cardinal::spatialShrunkenCentroids(CardinalObj, ...)
# new_CardinalObj <- add_ssc_annotation(CardinalObj, ssc_data, resolution ="r=2,k=8,s=32")
```
