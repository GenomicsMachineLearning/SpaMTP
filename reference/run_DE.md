# Runs EdgeR analysis for pooled data

Worker function for calculating differentially abundant metabolites per
pooling group. This function is used by by
[`FindAllDEMs()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindAllDEMs.md).

## Usage

``` r
run_DE(
  pooled_data,
  seurat_data,
  ident,
  output_dir,
  run_name,
  n,
  logFC_threshold,
  annotation.column,
  assay,
  return.individual = FALSE,
  verbose = TRUE
)
```

## Arguments

- pooled_data:

  A SingleCellExperiment object which contains the pooled
  pseudo-replicate data.

- seurat_data:

  A Seurat object containing the merged Xenium data being analysed (this
  is subset).

- ident:

  A character string defining the ident column to perform differential
  expression analysis against.

- output_dir:

  A character string defining the ident column to perform differential
  expression analysis against.

- run_name:

  A character string defining the title of this DE analysis (will be
  used when saving DEMs to .csv file).

- n:

  An integer that defines the number of pseudo-replicates per sample
  (default = 3).

- logFC_threshold:

  A numeric value indicating the logFC threshold to use for defining
  significant genes (default = 1.2).

- annotation.column:

  Character string defining the column where annotation information is
  stored in the assay metadata. This requires AnnotateSeuratMALDI() to
  be run where the default column to store annotations is
  "all_IsomerNames" (default = "None").

- assay:

  A character string defining the assay where the mz count data and
  annotations are stored (default = "Spatial").

- return.individual:

  Boolean value defining whether to return a list of individual edgeR
  objects for each designated ident. If FALSE, one merged edgeR object
  will be returned (default = FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A modified edgeR object which contains the relative pseudo-bulking
analysis outputs, including a DEMs data.frame with a list of
differential expressed m/z metabolites

## Examples

``` r
# pooled_obj <- run_pooling(SeuratObj, "sample", n = 3)
# run_DE(pooled_obj, SeuratObj, "sample", "~/Documents/DE_output/", "run_1", n = 3, logFC_threshold = 1.2, annotation.column = "all_IsomerNames", assay = "Spatial")
```
