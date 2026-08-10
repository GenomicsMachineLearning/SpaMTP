# Annotate one or more observed m/z values

Convenience wrapper that builds an index when one is not supplied. Reuse
a pre-built index for repeated or large searches.

## Usage

``` r
AnnotateMZ(
  observed_mz,
  db = NULL,
  index = NULL,
  polarity = c("positive", "negative", "neutral"),
  adducts = NULL,
  rules = NULL,
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

  Ion mode.

- adducts:

  Optional adduct subset.

- rules:

  Optional custom rule table.

- ppm:

  Mass tolerance in ppm.

- ms1_spectrum:

  Optional contextual spectrum.

- ...:

  Additional arguments passed to
  [`QueryMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/QueryMZAnnotationIndex.md).

## Value

A ranked candidate data frame.
