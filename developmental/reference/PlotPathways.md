# Plots the expression profile of a feature set corresponding to specified pathways onto a 2D scatter plot based on a dimensionality reduction technique.

This function has been adapted from fgsea package. For more detail about
function inputs please visit their
[documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).

## Usage

``` r
PlotPathways(
  pathways,
  object,
  title = NULL,
  assay = DefaultAssay(object),
  slot = "scale.data",
  reduction = NULL,
  colors = c("darkblue", "lightgrey", "darkred"),
  guide = "colourbar",
  ...
)
```

## Arguments

- pathways:

  Character vector of pathway names to plot. These pathways must be
  those stored in the `chempathway` object. To check call
  `unique(chempathway$pathwayName)`.

- object:

  SpaMTP Seurat object containing the dimensionality reduction to use
  for plotting.

- title:

  Optional title for the plot. If set to `NULL` the title will be the
  pathway name (default = NULL).

- assay:

  Character string defining the name of assay to use for data
  extraction. This `assay` should contain RAMP_IDs as feature names
  (default = `SeuratObject::DefaultAssay(object)`).

- slot:

  Character string defining the name of slot to extract data from
  (default = "scale.data").

- reduction:

  Character string defining which dimensionality reduction to use. If
  not specified, first searches for 'umap', then 'tsne', then 'pca'
  (default = NULL).

- colors:

  Vector of colors for the gradient (default = c("darkblue",
  "lightgrey", "darkred")).

- guide:

  Character string stating the type of legend to display (default =
  "colourbar").

- ...:

  Additional inputs taken by
  [`Seurat::FeaturePlot()`](https://satijalab.org/seurat/reference/FeaturePlot.html).
  Check the relative documentation for more infomation.

## Value

A ggplot object visualizing the score of each pathway across the
relative reduction.

## Examples

``` r
#PlotPathways(c("Glycolysis", "Acylcarnitine 3-Butenylcarnitine", "ABC transporters"), spamtp_obj, reduction = "umap"))
```
