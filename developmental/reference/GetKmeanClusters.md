# Perform K-means clustering on a specified reduction

This function runs K-means clustering on a specified reduction in a
SpaMTP Seurat object and adds the cluster assignments to the object
metadata.

## Usage

``` r
GetKmeanClusters(
  data,
  reduction = "SpatialPCA",
  cluster.name = "spatial_clusters",
  clusters = 8,
  iter.max = 10,
  nstart = 1,
  algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
  trace = FALSE,
  seed = 888
)
```

## Arguments

- data:

  A SpaMTP Seurat object containing the results from
  [`RunSpatialGraphPCA()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/RunSpatialGraphPCA.md).

- reduction:

  Character string stating the name of the reduction slot to use
  (default = "SpatialPCA").

- cluster.name:

  Character string of the name of the metadata column to store the
  cluster labels (default = "spatial_clusters").

- clusters:

  Integer defining the number of clusters to form (default = 8).

- iter.max:

  Integer defining the maximum number of iterations allowed (default =
  10).

- nstart:

  Integer stating the number of random sets to choose (default = 1).

- algorithm:

  Character string defining the K-means algorithm to use. One of
  `"Hartigan-Wong"`, `"Lloyd"`, `"Forgy"`, or `"MacQueen"` (default =
  "Hartigan-Wong").

- trace:

  Logical boolean indicating whether to produce tracing information on
  the progress of the algorithm (default = FALSE).

- seed:

  Integer of the random seed to use for reproducibility (default = 888).

## Value

A SpaMTP Seurat object with a new metadata column containing the K-means
cluster assignments.

## Examples

``` r
# seurat_object <- GetKmeanClusters(spamtp_obj, reduction = "SpatialPCA", centers = 8, cluster.name = "test_clusters")
# SpatialDimPlot(spamtp_obj, group.by = "test_clusters")
```
