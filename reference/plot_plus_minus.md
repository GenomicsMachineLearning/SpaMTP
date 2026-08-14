# Helper Function to generate merged counts within the plus minus range provided

Helper Function to generate merged counts within the plus minus range
provided

## Usage

``` r
plot_plus_minus(object, mz_list, plusminus)
```

## Arguments

- object:

  Seurat Spatial Metabolomic object containing mz values.

- mz_list:

  Vector of numeric m/z values to plot (e.g. c(mz-400.1578, mz-300.1)).

- plusminus:

  Numeric value defining the range/threshold either side of each mz peak
  provided.

## Value

List contatining a Seurat data object which contains the binned counts
in the `@metadata` slot, the column name under where the data is stored
and the title of this binned which will be plotted.

## Examples

``` r
### HELPER FUNCTION ###
```
