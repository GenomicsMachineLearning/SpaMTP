# Manually align an image (e.g. H&E, Immuno) to a SM SpaMTP dataset

Manually align an image (e.g. H&E, Immuno) to a SM SpaMTP dataset

## Usage

``` r
AddSMImage(
  image_path,
  SpaMTP,
  fov = "fov",
  grey.scale = 0.5,
  plot.greyscale = FALSE,
  seed = 123,
  n.spots = NULL,
  ...
)
```

## Arguments

- image_path:

  Character string defining the full path of the image to align.

- SpaMTP:

  SpaMTP Seurat object to align the image to.

- fov:

  Character string defining the image fov that contains the SM data
  coordinates (default = "fov").

- grey.scale:

  Numeric value defining the grey scale to use for generating a tissue
  mask from the provided image (default = 0.5).

- plot.greyscale:

  Boolean indicating whether to display the grey scale plot used to
  generate the tissue mask (default = FALSE).

- seed:

  Integer value defining the seed to use for calculating random fake
  gene values for aligning the image (default = 123).

- n.spots:

  Integer specifying the number of fake spots to generate for the
  aligned image. If NULL the number of spots will match that of the
  SpaMTP object provided (default = NULL).

- ...:

  Additional inputs used by the AlignSpatialOmics function. Please see
  documentation or call ?AlignSpatialOmics for more infomation.

## Value

A SpaMTP Seurat object with the provided image aligned to the SM data.
The image is stored in the `@image$slice1` slot.

## Examples

``` r
# AddSMImage(image_path = "../HnE_image.png", SpaMTP = SpaMTP_obj)
```
