# Add SMILES-derived structural features to a metabolite database

Add SMILES-derived structural features to a metabolite database

## Usage

``` r
AnnotateSMILESStructure(
  db,
  smiles_column = NULL,
  backend = c("auto", "native"),
  overwrite = FALSE,
  strict = FALSE,
  workers = getOption("SpaMTP.smiles_workers", 1L)
)
```

## Arguments

- db:

  A metabolite data frame.

- smiles_column:

  Column containing SMILES. When `NULL`, SpaMTP detects `iso_smiles`,
  `canonical_smiles`, or `smiles`.

- backend:

  Parser backend passed to
  [`DeconvolveSMILES()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/DeconvolveSMILES.md).

- overwrite:

  Replace existing feature values. By default only absent or missing
  values are filled.

- strict:

  Passed to
  [`DeconvolveSMILES()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/DeconvolveSMILES.md).

- workers:

  Parallel workers passed to
  [`DeconvolveSMILES()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/DeconvolveSMILES.md).

## Value

`db` with structure-derived columns appended or completed.

## Examples

``` r
utils::str(formals(AnnotateSMILESStructure))
#> Dotted pair list of 6
#>  $ db           : symbol 
#>  $ smiles_column: NULL
#>  $ backend      : language c("auto", "native")
#>  $ overwrite    : logi FALSE
#>  $ strict       : logi FALSE
#>  $ workers      : language getOption("SpaMTP.smiles_workers", 1L)
```
