# Calculates the mz range of the observed_df

Calculates the mz range of the observed_df

## Usage

``` r
calculate_bounds(input_df)
```

## Arguments

- input_df:

  DataFrame of the observed dataframe being annotated

## Value

A list containing the lower and upper mz range for the provided sample

## Examples

``` r
# mz_df <- SeuratObject`[["Spatial"]][["mz"]]`
# mz_df$row_id <- seq(1, length(mz_df$mz))

# mass_range <- calculate_bounds(mz_df)
# lower_bound <- mass_range$lower_bound
# upper_bound <- mass_range$upper_bound
```
