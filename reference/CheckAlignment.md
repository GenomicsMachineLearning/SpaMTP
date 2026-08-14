# Check multi-modal coordinate alignment

Checks the alignment of two spatial datasets by plotting their relative
coordinates on the same graph

## Usage

``` r
CheckAlignment(
  SM.data,
  ST.data,
  image.res = NULL,
  names = c("SM", "ST"),
  cols = NULL,
  image.slice = "slice1",
  size = 0.5
)
```

## Arguments

- SM.data:

  SpaMTP Seurat object containing SM data

- ST.data:

  SpaMTP Seurat object containing ST data

- image.res:

  Character string defining the Visium image resolution to use. This is
  required for the correct scale.factor to be applied (default = NULL).

- names:

  Vector of 2 character strings used to define each dataset being
  plotted (default = c("SM", "ST")).

- cols:

  Vector of 2 colors used for plotting (default = NULL).

- image.slice:

  Character string matching the image slice name within the ST SpaMTP
  Seurat object (default = "slice1").

- size:

  Numeric value indicating the point size to plot (default = 0.5).

## Value

A 2D scatter plot showing the relative spatial locations of the ST and
SM data points

## Examples

``` r
# CheckAlignment(SM.data, ST.data)
```
