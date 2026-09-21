# Metaspace R Client Constructor

Creates a client object with methods to interact with the METASPACE
GraphQL API.

## Usage

``` r
metaspace_client(host = "https://metaspace2020.org", api_key = NULL)
```

## Arguments

- host:

  Character string specifying the METASPACE host URL (default =
  "https://metaspace2020.org").

- api_key:

  Optional. Your METASPACE API key for accessing private datasets. If
  accessing a public dataset this can be left as `NULL` (default =
  NULL).

## Value

A list containing functions for querying METASPACE data.

## Examples

``` r
utils::str(formals(metaspace_client))
#> Dotted pair list of 2
#>  $ host   : chr "https://metaspace2020.org"
#>  $ api_key: NULL
# ms <- metaspace_client()
```
