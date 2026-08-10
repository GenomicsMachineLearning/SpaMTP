# Checks if the complete adduct is in the data base, else returns a truncated adduct

Checks if the complete adduct is in the data base, else returns a
truncated adduct

## Usage

``` r
check_and_truncate_adduct_vector(adduct, db, verbose = TRUE)
```

## Arguments

- adduct:

  Character string defining the adduct to be checked.

- db:

  DataFrame of the current reference database.

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = TRUE).

## Value

Character string of the complete or truncated adduct

## Examples

``` r
# HMDB_db <- load("data/HMDB_1_names.rds")
# check_and_truncate_adduct_vector(c("M+H"), HMDB_db)
```
