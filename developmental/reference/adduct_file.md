# Common adduct constants

A small table of legacy adduct definitions retained for compatibility
with earlier SpaMTP workflows. New annotation code should use
[`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md),
which includes charge, stoichiometry, ion-mode, and chemical-validity
fields.

## Usage

``` r
adduct_file
```

## Format

A data frame with 47 rows and 6 variables:

- adduct_name:

  Adduct notation.

- ion.mass:

  Legacy ion-mass expression.

- charge:

  Ion charge.

- mult:

  Analyte stoichiometric multiplier.

- add_mass:

  Exact mass shift.

- pol:

  Ion polarity.

## Value

A data frame of legacy adduct definitions.

## See also

[`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/AdductRules.md),
[`MALDIMatrixRules()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/MALDIMatrixRules.md)
