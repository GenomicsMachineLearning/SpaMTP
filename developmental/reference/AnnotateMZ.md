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

- ...:

  Additional arguments passed to
  [`QueryMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/QueryMZAnnotationIndex.md).

## Value

A ranked candidate data frame.
