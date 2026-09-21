# Runs multilevel Monte-Carlo variant for performing gene sets co-regulation analysis using the RAMP_DB metabolite/gene database.

This function is adapted from the
[fgsea::geseca](https://github.com/alserglab/fgsea/blob/master/R/geseca-multilevel.R)
package to identify significantly expressed RAMP_DB pathways based on an
expression/feature embedding matrix.

## Usage

``` r
RunRAMPgeseca(
  E,
  minSize = 1,
  maxSize = nrow(E) - 1,
  center = TRUE,
  scale = FALSE,
  sampleSize = 101,
  eps = 1e-50,
  nproc = 0,
  BPPARAM = NULL,
  nPermSimple = 1000,
  database = NULL,
  database_version = "latest",
  database_source = c("auto", "bundled", "local"),
  database_local_dir = NULL
)
```

## Arguments

- E:

  expression matrix, rows corresponds to RAMP_IDs, columns corresponds
  to cell barcodes.

- minSize:

  Minimal size of a gene set to test. All pathways below the threshold
  are excluded (default = 1).

- maxSize:

  Maximal size of a gene set to test. All pathways above the threshold
  are excluded (default = `nrow(E) - 1`).

- center:

  a logical value indicating whether the gene expression should be
  centered to have zero mean before the analysis takes place (default =
  TRUE).

- scale:

  a logical value indicating whether the gene expression should be
  scaled to have unit variance before the analysis takes place (default
  = FALSE).

- sampleSize:

  sample size for conditional sampling (default = 101).

- eps:

  This parameter sets the boundary for calculating P-values (default =
  1e-50).

- nproc:

  If not equal to zero sets BPPARAM to use nproc workers (default = 0).

- BPPARAM:

  Parallelization parameter used in bplapply (default = NULL).

- nPermSimple:

  Number of permutations in the simple geseca implementation for
  preliminary estimation of P-values (default = 1000).

- database:

  Optional named list of database resources, normally created by
  [`LoadSpaMTPDatabase()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/LoadSpaMTPDatabase.md).

- database_version:

  Database snapshot version used for pathway lookup.

- database_source:

  Database source; see
  [`LoadSpaMTPDatabase()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/LoadSpaMTPDatabase.md).

- database_local_dir:

  Optional local RDS resource directory.

## Value

A table with GESECA results. Each row corresponds to a tested RAMP_DB
pathway.

## Examples

``` r
utils::str(formals(RunRAMPgeseca))
#> Dotted pair list of 14
#>  $ E                 : symbol 
#>  $ minSize           : num 1
#>  $ maxSize           : language nrow(E) - 1
#>  $ center            : logi TRUE
#>  $ scale             : logi FALSE
#>  $ sampleSize        : num 101
#>  $ eps               : num 1e-50
#>  $ nproc             : num 0
#>  $ BPPARAM           : NULL
#>  $ nPermSimple       : num 1000
#>  $ database          : NULL
#>  $ database_version  : chr "latest"
#>  $ database_source   : language c("auto", "bundled", "local")
#>  $ database_local_dir: NULL
# E <- SpaMTP@reductions$pca.rev@feature.loadings
# sig_pathways <- RunRAMPgeseca(E, minSize=15, maxSize=500)
```
