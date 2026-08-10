# Find top features and metabolites that are strongly correlated with a given feature

This function returns a list of features ranked by highest Pearson
correlation score. The specified feature to corrleate against can be
either a m/z value, gene or an ident (i.e. cluster). For multi-omic
data, both a metabolic and transcriptomic assay can be specified to
calculate correlation of both metabolites and genes.

## Usage

``` r
FindCorrelatedFeatures(
  data,
  mz = NULL,
  gene = NULL,
  ident = NULL,
  SM.assay = "SPM",
  ST.assay = NULL,
  SM.slot = "counts",
  ST.slot = "counts",
  nfeatures = 10
)
```

## Arguments

- data:

  SpaMTP Seurat class object containing both Spatial Transcriptomic and
  Metabolic data assays.

- mz:

  Numeric string specifying the m/z to find correlated features for. One
  of `mz`, `gene` or `ident` must be provided, alternatives must be
  `NULL` (default = NULL).

- gene:

  Character string specifying the gene to find correlated features for.
  One of `mz`, `gene` or `ident` must be provided, alternatives must be
  `NULL` (default = NULL).

- ident:

  Character string defining the ident column in the data object's
  `@meta.data` slot to find correlated features for. One of `mz`, `gene`
  or `ident` must be provided, alternatives must be `NULL` (default =
  NULL).

- SM.assay:

  Character string specifying the name of the assay containing the
  spatial metabolomics (SM) data (default = "SPM").

- ST.assay:

  Character string specifying the name of the assay containing the
  spatial transcriptomics (ST) data. If NULL then only metabolites will
  be used (Default = NULL).

- SM.slot:

  Character string specifying the slot of the SM assay to use (default =
  "counts").

- ST.slot:

  Character string specifying the slot of the ST assay to use (default =
  "counts").

- nfeatures:

  Integer specifying the number of top correlated features to return
  (default = 10).

## Value

A data frame containing the top correlated features with columns for the
feature names and their correlation values.

## Examples

``` r
# result <- FindCorrelatedFeatures(data = SpaMTP, gene = "GeneX", nfeatures = 5)
```
