# Filters reference Database to only select natural elements

Filters reference Database to only select natural elements

## Usage

``` r
formula_filter(df, elements = NULL)
```

## Arguments

- df:

  DataFrame of the reference database.

- elements:

  Vector of character strings of elements to include (default = c("H",
  "C", "N", "O", "S", "Cl", "Br", "F", "Na", "P", "I", "Si")).

## Value

A refined DataFrame which only includes annotations containing the
specified elements

## Examples

``` r
## Get filtered DB by adduct
# db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos")

## Refine to only include natural elements
# db_2 <- formula_filter(db_1)
```
