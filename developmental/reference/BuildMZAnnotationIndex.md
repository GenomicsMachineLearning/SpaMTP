# Build a reusable indexed metabolite annotation search space

Expected adduct m/z values are generated once, sorted, and queried by
binary interval search. This is the one-dimensional equivalent of an
interval tree and avoids scanning every database row for every observed
peak.

## Usage

``` r
BuildMZAnnotationIndex(
  db,
  polarity = c("positive", "negative", "neutral"),
  adducts = NULL,
  rules = NULL,
  collapse_isomers = TRUE
)
```

## Arguments

- db:

  Metabolite database. Both the legacy SpaMTP five-column/wide database
  and the RaMP `chem_props` schema are supported.

- polarity:

  `"positive"`, `"negative"`, or `"neutral"`.

- adducts:

  Optional character vector of adduct names or bracketed notations.

- rules:

  Optional custom rule table. See
  [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md).

- collapse_isomers:

  Collapse records sharing formula, exact mass, and proton bound before
  indexing.

## Value

An object of class `spamtp_mz_index`.
