# Annotates FMP10 matrix data

Adds metabolite annotations to respective m/z values generated using an
FMP10 matrix. This is done based on a curated FMP10 matrix database.

## Usage

``` r
AddFMP10Annotations(
  obj,
  only.fmp.adduct = FALSE,
  add.custom.annotation = NULL,
  assay = "Spatial",
  return.only.annotated = FALSE,
  mass.threshold = 0.05,
  annotation.column = "all_IsomerNames",
  database = NULL,
  database_version = "latest",
  database_source = c("auto", "bundled", "local"),
  database_local_dir = NULL
)
```

## Arguments

- obj:

  SpaMTP Seurat object containing m/z intensity values for annotation.
  Data should be generated with a FMP10 matrix.

- only.fmp.adduct:

  Boolean indicating if only metabolites with FMP10+ adducts (`+FMP10`,
  +`2FMP10`, etc.) should be assigned to m/z values (default = FALSE).

- add.custom.annotation:

  data.frame containing addition FMP10 metabolite annotations that are
  not in the current FMP10 database. Note: this data.frame must contain
  these column c("mass", "annotation", "Adduct", "Formula", "Isomers",
  "Isomers_IDs"). If set to NULL, only reference FMP10 database will be
  used (default = NULL).

- assay:

  Character string defining the Seurat object assay to store the
  respective annotations in the feature meta.data dataframe (default =
  "Spatial").

- return.only.annotated:

  Boolean defining whether to return a SpaMTP Seurat object containing
  only successfully annotated m/z values (default = FALSE).

- mass.threshold:

  Numeric value defining the acceptable threshold (plus-minus) between
  the custom annotations and the actual m/z values contained within the
  SpaMTP object (default = 0.05).

- annotation.column:

  Character string defining the feature meta.data column name that will
  contain the assigned annotations (default = "all_IsomerNames").

- database:

  Optional named list of database resources, normally created by
  [`LoadSpaMTPDatabase()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/LoadSpaMTPDatabase.md).

- database_version:

  Database snapshot version used for annotation lookup.

- database_source:

  Database source; see
  [`LoadSpaMTPDatabase()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/LoadSpaMTPDatabase.md).

- database_local_dir:

  Optional local RDS resource directory.

## Value

SpaMTP Seurat object containing the relative metabolite annotations
stored in the feature metadata dataframe.

## Examples

``` r
utils::str(formals(AddFMP10Annotations))
#> Dotted pair list of 11
#>  $ obj                  : symbol 
#>  $ only.fmp.adduct      : logi FALSE
#>  $ add.custom.annotation: NULL
#>  $ assay                : chr "Spatial"
#>  $ return.only.annotated: logi FALSE
#>  $ mass.threshold       : num 0.05
#>  $ annotation.column    : chr "all_IsomerNames"
#>  $ database             : NULL
#>  $ database_version     : chr "latest"
#>  $ database_source      : language c("auto", "bundled", "local")
#>  $ database_local_dir   : NULL
# AddFMP10Annotations(spamtp, only.fmp.adduct = FALSE)
```
