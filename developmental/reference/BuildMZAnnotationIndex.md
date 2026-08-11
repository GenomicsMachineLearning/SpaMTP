# Build a reusable indexed metabolite annotation search space

Expected adduct m/z values are generated once, sorted, and queried by
binary interval search. This is the one-dimensional equivalent of an
interval tree and avoids scanning every database row for every observed
peak.

## Usage

``` r
BuildMZAnnotationIndex(
  db,
  polarity = NULL,
  adducts = NULL,
  rules = NULL,
  maldi_matrix = NULL,
  collapse_isomers = TRUE
)
```

## Arguments

- db:

  Metabolite database. Both the legacy SpaMTP five-column/wide database
  and the RaMP `chem_props` schema are supported.

- polarity:

  `"positive"`, `"negative"`, or `"neutral"`. When `NULL`, use the
  matrix-profile default, or positive mode when no profile is given.

- adducts:

  Optional character vector of adduct names or bracketed notations. When
  `NULL`, use the complete automatically selected rule space; this is a
  filter, not a compulsory input.

- rules:

  Optional custom rule table. See
  [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md).

- maldi_matrix:

  Optional MALDI matrix or derivatization reagent profile. When supplied
  and `rules` is `NULL`,
  [`MALDIMatrixRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/MALDIMatrixRules.md)
  automatically selects standard and validated matrix-specific rules.
  `adducts` remains an optional filter on that selected rule space.

- collapse_isomers:

  Collapse records sharing formula, exact mass, and proton bound before
  indexing.

## Value

An object of class `spamtp_mz_index`.
