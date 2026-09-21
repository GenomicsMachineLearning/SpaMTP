# Annotate one or more observed m/z values

Convenience wrapper that builds an index when one is not supplied. Reuse
a pre-built index for repeated or large searches.

## Usage

``` r
AnnotateMZ(
  observed_mz,
  db = NULL,
  index = NULL,
  polarity = NULL,
  adducts = NULL,
  rules = NULL,
  maldi_matrix = NULL,
  ppm = 5,
  ms1_spectrum = NULL,
  database_version = "latest",
  database_source = c("auto", "bundled", "local"),
  database_local_dir = NULL,
  infer_structure = c("auto", "never", "always"),
  structure_backend = c("auto", "native"),
  structure_workers = getOption("SpaMTP.smiles_workers", 1L),
  min_structure_score = 0.05,
  ...
)
```

## Arguments

- observed_mz:

  Numeric vector of observed m/z values.

- db:

  Metabolite database used when `index` is `NULL`.

- index:

  Optional pre-built `spamtp_mz_index`.

- polarity:

  Ion mode. When `NULL`, infer it from `index` or the matrix profile,
  falling back to positive mode.

- adducts:

  Optional adduct subset. When `NULL`, retain the complete rule space
  selected by the matrix profile or polarity.

- rules:

  Optional custom rule table.

- maldi_matrix:

  Optional MALDI matrix/reagent profile used to select rules
  automatically when `rules` is `NULL`.

- ppm:

  Mass tolerance in ppm.

- ms1_spectrum:

  Optional contextual spectrum.

- database_version:

  Database snapshot version used when `db = NULL`.

- database_source:

  Database source used when `db = NULL`; see
  [`LoadSpaMTPDatabase()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/LoadSpaMTPDatabase.md).

- database_local_dir:

  Optional local RDS resource directory.

- infer_structure, structure_backend, structure_workers,
  min_structure_score:

  Structure-aware rule-selection arguments passed to
  [`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/BuildMZAnnotationIndex.md).

- ...:

  Additional arguments passed to
  [`QueryMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/QueryMZAnnotationIndex.md).

## Value

A ranked candidate data frame.

## Examples

``` r
utils::str(formals(AnnotateMZ))
#> Dotted pair list of 17
#>  $ observed_mz        : symbol 
#>  $ db                 : NULL
#>  $ index              : NULL
#>  $ polarity           : NULL
#>  $ adducts            : NULL
#>  $ rules              : NULL
#>  $ maldi_matrix       : NULL
#>  $ ppm                : num 5
#>  $ ms1_spectrum       : NULL
#>  $ database_version   : chr "latest"
#>  $ database_source    : language c("auto", "bundled", "local")
#>  $ database_local_dir : NULL
#>  $ infer_structure    : language c("auto", "never", "always")
#>  $ structure_backend  : language c("auto", "native")
#>  $ structure_workers  : language getOption("SpaMTP.smiles_workers", 1L)
#>  $ min_structure_score: num 0.05
#>  $ ...                : symbol 
```
