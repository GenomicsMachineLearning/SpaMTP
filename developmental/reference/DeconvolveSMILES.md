# Decompose explicit functional groups in SMILES strings

A lightweight, dependency-free SMILES graph parser identifies common
metabolite functional groups and derives interpretable positive-ion,
negative-ion, alkali-binding, and neutral-mass priors. The parser is
intended for candidate prioritisation rather than pKa, proton-affinity,
or quantum chemistry calculations. Stereochemistry is accepted but does
not change the counts.

## Usage

``` r
DeconvolveSMILES(
  smiles,
  backend = c("auto", "native"),
  strict = FALSE,
  workers = getOption("SpaMTP.smiles_workers", 1L)
)
```

## Arguments

- smiles:

  Character vector of SMILES strings.

- backend:

  Currently `"auto"` and `"native"` both use SpaMTP's native parser. The
  argument reserves a stable extension point for optional chemistry
  toolkits.

- strict:

  Stop when any SMILES cannot be parsed. When `FALSE`, invalid rows are
  returned with `structure_valid = FALSE` and missing features.

- workers:

  Number of forked workers used for unique SMILES on platforms
  supporting
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).
  The default is controlled by `options(SpaMTP.smiles_workers = 1)`.

## Value

A data frame with functional-group counts, ion-mode scores, and a
human-readable `structure_evidence` field.

## Examples

``` r
utils::str(formals(DeconvolveSMILES))
#> Dotted pair list of 4
#>  $ smiles : symbol 
#>  $ backend: language c("auto", "native")
#>  $ strict : logi FALSE
#>  $ workers: language getOption("SpaMTP.smiles_workers", 1L)
```
