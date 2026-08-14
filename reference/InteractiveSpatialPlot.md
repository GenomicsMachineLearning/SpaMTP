# Interactive Spatial Plot for visualising different m/z bin sizes

This function implements a shinyApp for interactive viusalisation of the
samples mean mass intensity spectra. Users can adjust the resolution/bin
size and observe changes in spatial expression respectively.

## Usage

``` r
InteractiveSpatialPlot(
  obj,
  assay = "Spatial",
  slot = "counts",
  image = "slice1"
)
```

## Arguments

- obj:

  SpaMTP object containing the intensity data to plot.

- assay:

  Character defining the Seurat object assay to extract the intensity
  data from (default = "Spatial").

- slot:

  Character defining the Seurat object assay slot to extract the
  intensity data from (default = "counts").

- image:

  Character defining the the name of the image to use for tissue
  coordinates and spatial plotting (default = "slice1").

## Examples

``` r
#InteractiveSpatialPlot(spamtp)
```
