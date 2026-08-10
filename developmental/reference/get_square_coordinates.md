# This function computes the coordinates of a square's four corners based on a given center point and width.

This function computes the coordinates of a square's four corners based
on a given center point and width.

## Usage

``` r
get_square_coordinates(center_x, center_y, width, name)
```

## Arguments

- center_x:

  Numeric value defining the x-coordinate of the square's center.

- center_y:

  Numeric value defining the y-coordinate of the square's center.

- width:

  Numeric value indicating the width of the square.

- name:

  Character string specifying the label assigned to the square for
  identification.

## Value

A data frame with columns `Selection`, `X`, and `Y` representing the
square's name and corner coordinates.

## Examples

``` r
get_square_coordinates(center_x = 5, center_y = 5, width = 4, name = "MySquare")
#>   Selection X Y
#> 1  MySquare 3 3
#> 2  MySquare 7 3
#> 3  MySquare 7 7
#> 4  MySquare 3 7
#> 5  MySquare 3 3
```
