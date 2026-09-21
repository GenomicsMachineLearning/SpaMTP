# Helper function for converting Seurat Class ggplots from spot to pixel layout

Helper function for converting Seurat Class ggplots from spot to pixel
layout

## Usage

``` r
pixelPlot(plot)
```

## Arguments

- plot:

  ggplot object contating the doplot to be converted into pixel layout

## Value

plot where spots are in pixel layout rather then spot

## Examples

``` r
utils::str(formals(pixelPlot))
#> Dotted pair list of 1
#>  $ plot: symbol 
#pixelPlot(SpatialFeaturPlot(SpaMTP.obj, features = "nFeature_Spatial"))
```
