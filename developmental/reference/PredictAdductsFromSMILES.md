# Predict a structure-aware adduct search space from SMILES

This function applies the same functional-group and ion-mode logic used
by
[`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/BuildMZAnnotationIndex.md)
before any observed m/z is supplied. It is useful for auditing which
protonation, deprotonation, alkali-binding, or reactive matrix
hypotheses SpaMTP will retain for a structure.

## Usage

``` r
PredictAdductsFromSMILES(
  smiles,
  polarity = NULL,
  maldi_matrix = NULL,
  rules = NULL,
  min_structure_score = 0.05,
  backend = c("auto", "native"),
  workers = getOption("SpaMTP.smiles_workers", 1L)
)
```

## Arguments

- smiles:

  Character vector of SMILES strings.

- polarity:

  `"positive"`, `"negative"`, or `"neutral"`. When `NULL`, use the
  matrix-profile default or positive mode.

- maldi_matrix:

  Optional MALDI matrix/reagent profile.

- rules:

  Optional custom rule table. It cannot be combined with `maldi_matrix`.

- min_structure_score:

  Minimum rule-specific structural prior marked as retained.

- backend, workers:

  Passed to
  [`DeconvolveSMILES()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/DeconvolveSMILES.md).

## Value

A ranked data frame with one row per structure/rule combination.

## Examples

``` r
utils::str(formals(PredictAdductsFromSMILES))
#> Dotted pair list of 7
#>  $ smiles             : symbol 
#>  $ polarity           : NULL
#>  $ maldi_matrix       : NULL
#>  $ rules              : NULL
#>  $ min_structure_score: num 0.05
#>  $ backend            : language c("auto", "native")
#>  $ workers            : language getOption("SpaMTP.smiles_workers", 1L)
```
