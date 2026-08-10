# Heatmap of Differentially Expressed Metabolites

Generates a heatmap of DEMs generated from edgeR analysis run using
[`FindAllDEMs()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/FindAllDEMs.md).
This function uses `pheatmap` to plot data.

## Usage

``` r
DEMsHeatmap(
  edgeR_output,
  n = 5,
  only.pos = FALSE,
  FDR.threshold = 0.05,
  logfc.threshold = 0.5,
  order.by = "FDR",
  scale = "row",
  color = (grDevices::colorRampPalette(c("navy", "white", "red")))(50),
  cluster_cols = F,
  cluster_rows = T,
  fontsize_row = 15,
  fontsize_col = 15,
  cutree_cols = 9,
  silent = TRUE,
  plot_annotations_column = NULL,
  save_to_path = NULL,
  plot.save.width = 20,
  plot.save.height = 20,
  nlabels.to.show = NULL,
  annotation_colors = NULL
)
```

## Arguments

- edgeR_output:

  A list containing outputs from edgeR analysis (from FindAllDEMs()).
  This includes pseudo-bulked counts and DEMs.

- n:

  A numeric integer that defines the number of UP and DOWN regulated
  peaks to plot (default = 25).

- only.pos:

  Boolean indicating if only positive markers should be returned
  (default = FALSE).

- FDR.threshold:

  Numeric value that defines the FDR threshold to use for defining most
  significant results (default = 0.05).

- logfc.threshold:

  Numeric value that defines the logFC threshold to use for filtering
  significant results (default = 0.5).

- order.by:

  Character string defining which parameter to order markers by, options
  are either 'FDR' or 'logFC' (default = "FDR").

- scale:

  A character string indicating if the values should be centered and
  scaled in either the row direction or the column direction, or none.
  Corresponding values are "row", "column" and "none"

- color:

  A vector of colors used in heatmap (default =
  grDevices::colorRampPalette(c("navy", "white", "red"))(50)).

- cluster_cols:

  Boolean value determining if columns should be clustered or hclust
  object (default = F).

- cluster_rows:

  Boolean value determining if rows should be clustered or hclust object
  (default = T).

- fontsize_row:

  A numeric value defining the fontsize of rownames (default = 15).

- fontsize_col:

  A numeric value defining the fontsize of colnames (default = 15).

- cutree_cols:

  A numeric value defining the number of clusters the columns are
  divided into, based on the hierarchical clustering(using cutree), if
  cols are not clustered, the argument is ignored (default = 9).

- silent:

  Boolean value indicating if the plot should not be draw (default =
  TRUE).

- plot_annotations_column:

  Character string indicating the column name that contains the
  metabolite annotations to plot. Annotations = TRUE must be used in
  FindAllDEMs() for edgeR output to include annotations. If
  plot_annotations_column = NULL, m/z vaues will be plotted (default =
  NULL).

- save_to_path:

  Character string defining the full filepath and name of the plot to be
  saved as.

- plot.save.width:

  Integer value representing the width of the saved pdf plot (default =
  20).

- plot.save.height:

  Integer value representing the height of the saved pdf plot (default =
  20).

- nlabels.to.show:

  Numeric value defining the number of annotations to show per m/z
  (default = NULL).

- annotation_colors:

  List for specifying annotation_row and annotation_col track colors
  manually. Check pheatmap R-Package documentation for details. If set
  to 'NA', default coloring will be used (default = NA).

## Value

A heatmap plot of significantly differentially expressed metabolites
defined in the edgeR ouput object.

## Examples

``` r
# DEMs <- FindAllDEMs(SeuratObj, "sample")

# DEMsHeatmap(DEMs)
```
