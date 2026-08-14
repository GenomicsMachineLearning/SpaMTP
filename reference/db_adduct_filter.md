# Filters the provided metabolomic database by polarity and adducts

Filters the provided metabolomic database by polarity and adducts

## Usage

``` r
db_adduct_filter(db, adduct, polarity = "neg", verbose = TRUE)
```

## Arguments

- db:

  DataFrame of the current reference database.

- adduct:

  Character string defining the adduct to be checked.

- polarity:

  Character string defining the polarity of the adducts (default =
  "neg").

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

A filtered reference metabolomic database DataFrame

## Examples

``` r
# HMDB_db <- load("data/HMDB_1_names.rds")
# db_adduct_filter(HMDB_db,c("M+H"), polarity = "pos")
```
