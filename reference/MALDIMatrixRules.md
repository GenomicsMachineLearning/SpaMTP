# Select annotation rules from the MALDI matrix/reagent profile

Conventional matrices select ordinary ion rules and, where validated,
matrix-metabolite adducts. Reactive matrices additionally select
covalent product rules. Profiles without a verified universal net
transformation are retained in the registry but do not generate
speculative mass shifts.

## Usage

``` r
MALDIMatrixRules(
  maldi_matrix,
  polarity = NULL,
  include_standard = TRUE,
  include_matrix_products = TRUE
)
```

## Arguments

- maldi_matrix:

  Matrix or derivatization reagent name; aliases such as `"2,5-DHB"`,
  `"HCCA"`, `"9-AA"`, and `"FMP-10"` are accepted.

- polarity:

  `"positive"`, `"negative"`, or `"neutral"`. If `NULL`, use the profile
  default; profiles supporting both modes require an explicit choice.

- include_standard:

  Include the validated ordinary adduct rules.

- include_matrix_products:

  Include validated matrix-adduct or reactive product rules.

## Value

A validated adduct/reaction rule data frame accepted by
[`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BuildMZAnnotationIndex.md).
