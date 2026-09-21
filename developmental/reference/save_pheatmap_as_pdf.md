# Saves a DEMsHeatmap as a PDF

Saves a DEMsHeatmap as a PDF

## Usage

``` r
save_pheatmap_as_pdf(pheatmap, filename, width = 20, height = 20)
```

## Arguments

- pheatmap:

  A pheatmap plot object that is being saved.

- filename:

  Character string defining the full filepath and name of the plot to be
  saved as.

- width:

  Integer value representing the width of the saved pdf plot (default =
  20).

- height:

  Integer value representing the height of the saved pdf plot (default =
  20).

## Value

The generated PDF path, invisibly.

## Examples

``` r
utils::str(formals(save_pheatmap_as_pdf))
#> Dotted pair list of 4
#>  $ pheatmap: symbol 
#>  $ filename: symbol 
#>  $ width   : num 20
#>  $ height  : num 20
# save_pheatmap_as_pdf(pheatmap, filename = "/Documents/plots/pheatmap1")
```
