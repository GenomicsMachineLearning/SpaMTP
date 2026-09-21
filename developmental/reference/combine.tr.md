# Combines rigid tranformation matrices

Combines rigid tranformation matrices in the following order:
translation of points to origin (0, 0) -\> reflection of points -\>
rotation by alpha degrees and translation of points to new center

## Usage

``` r
combine.tr(center.cur, center.new, alpha, mirror.x = FALSE, mirror.y = FALSE)
```

## Arguments

- center.cur:

  (x, y) image pixel coordinates specifying the current center of the
  tissue (stored in slot "tools" as "centers")

- center.new:

  (x, y) image pixel coordinates specifying the new center (image
  center)

- alpha:

  Rotation angle

- mirror.x, mirror.y:

  Logical values indicating reflection across the x or y axis before
  rotation.

## Value

A 3-by-3 homogeneous affine-transformation matrix.

## Examples

``` r
utils::str(formals(combine.tr))
#> Dotted pair list of 5
#>  $ center.cur: symbol 
#>  $ center.new: symbol 
#>  $ alpha     : symbol 
#>  $ mirror.x  : logi FALSE
#>  $ mirror.y  : logi FALSE
transformation <- combine.tr(
  center.cur = c(0, 0),
  center.new = c(10, 20),
  alpha = 90
)
stopifnot(identical(dim(transformation), c(3L, 3L)))
```
