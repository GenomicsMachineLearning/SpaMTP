# Finds differentially expressed m/z values/metabolites between all comparison groups.

Finds differentially expressed m/z values/metabolites between all
comparison groups.

## Usage

``` r
FindAllDEMs(
  data,
  ident,
  n = 3,
  logFC_threshold = 1.2,
  DE_output_dir = NULL,
  run_name = "FindAllDEMs",
  annotation.column = NULL,
  assay = "Spatial",
  slot = "counts",
  return.individual = FALSE,
  verbose = TRUE,
  seed = 1234
)
```

## Arguments

- data:

  A Seurat object containing mz values for differential expression
  analysis.

- ident:

  A character string defining the metadata column or groups to compare
  mz values between.

- n:

  An integer that defines the number of pseudo-replicates (pools) per
  sample (default = 3).

- logFC_threshold:

  A numeric value indicating the logFC threshold to use for defining
  significant genes (default = 1.2).

- DE_output_dir:

  A character string defining the directory path for all output files to
  be stored. This path must a new directory. Else, set to NULL as
  default.

- run_name:

  A character string defining the title of this DE analysis that will be
  used when saving DEMs to .csv file (default = 'FindAllDEMs').

- annotation.column:

  Character string defining the column where annotation information is
  stored in the assay metadata. This requires AnnotateSeuratMALDI() to
  be run where the default column to store annotations is
  "all_IsomerNames" (default = "None").

- assay:

  A character string defining the assay where the mz count data and
  annotations are stored (default = "Spatial").

- slot:

  Character string defining the assay storage slot to pull the relative
  mz intensity values from. Note: EdgeR requires raw counts, all values
  must be positive (default = "counts").

- return.individual:

  Boolean value defining whether to return a list of individual edgeR
  objects for each designated ident. If FALSE, one merged edgeR object
  will be returned (default = FALSE).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

- seed:

  Numeric value used to set the seed for reproducible randomisation
  (default = 1234).

## Value

Returns an list() contains the EdgeR DE results. Pseudo-bulk counts are
stored in \$counts and DEMs are in \$DEMs.

## Examples

``` r
# FindAllDEMs(SeuratObj, "sample",DE_output_dir = "~/Documents/DE_output/", annotations = TRUE)
```
