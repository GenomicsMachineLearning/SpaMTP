# Construct a k-Nearest Neighbour Graph from Spatial Coordinates.

Builds a sparse adjacency matrix representing the k-nearest neighbour
relationships between spots or cells, using Euclidean distance between
spatial coordinates.

## Usage

``` r
kneighbors_graph(
  location,
  n_neighbors,
  platform,
  include_self = FALSE,
  alg = "kd_tree"
)
```

## Arguments

- location:

  A data frame or matrix with columns `x` and `y` representing spatial
  coordinates, and rownames match barcode/cell names.

- n_neighbors:

  Integer value defining the number of neighbours to consider.

- platform:

  Character string being either `"Visium"` or `"ST"`, specifying the
  platform of the dataset being analysed. Determines how distance
  thresholds and adjacency are calculated whereby "Visium" data is
  handled slightly differently to compensate for the hexagon shape of
  the spots.

- include_self:

  Boolean logical value stating whether to include self-connections
  (default = FALSE).

- alg:

  Character string specifying the algorithm to use for nearest neighbour
  search. Passed to
  [`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html) (default
  = "kd_tree").

## Value

A sparse adjacency matrix (`dgCMatrix`) representing the k-NN graph.

## Examples

``` r
### HELPER FUNCTION
```
