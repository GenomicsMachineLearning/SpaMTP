# Return the validated SpaMTP adduct rule table

The rules use neutral monoisotopic mass and explicit ion masses. Each
row stores the molecular multiplier, total numerator mass shift, charge,
number of protons removed from the analyte, and the charge contributed
by added species. This allows rules to be checked before candidate
generation.

## Usage

``` r
AdductRules(
  polarity = c("both", "positive", "negative", "neutral"),
  include_complex = TRUE
)
```

## Arguments

- polarity:

  One of `"both"`, `"positive"`, `"negative"`, or `"neutral"`.

- include_complex:

  Include multimers, solvent clusters, metal-exchange, and multiply
  charged ions.

## Value

A data frame containing validated adduct rules.
