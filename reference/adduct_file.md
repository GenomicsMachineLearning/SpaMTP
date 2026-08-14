# adduct_file: A dataframe containing possible adducts used for pathway analysis

This object contains a collection of adducts with their relative
ion.mass, change and polarity

## Usage

``` r
adduct_file
```

## Format

### A data frame with 47 rows and 6 variables:

- adduct_name:

  Chemical formula of the adduct (character)

- ion.mass:

  Mass of the molecular ion formula (character)

- charge:

  Relative charge of the adduct (integer)

- mult:

  Multiplication factor based on the original molecule mass (double)

- add_mass:

  The true mass of the ion being added or reduced (double)

- pol:

  Polarity of the ion (reduced = negative, added = positive) (character)
