# Interactive app for SM and ST coordinate alignment

Shiny app allowing for manual alignment of SM pixel coordinates to the
same coordinate system of the provided ST data.

## Usage

``` r
AlignSpatialOmics(
  sm.data,
  st.data,
  msi.pixel.multiplier = 20,
  image.res = "lowres",
  continous_cols = NULL,
  catagorical_cols = NULL,
  fov = "fov",
  image.slice = "slice1",
  shiny.host = "0.0.0.0",
  shiny.port = 4698,
  verbose = FALSE
)
```

## Arguments

- sm.data:

  SpaMTP Seurat Object containing SM data

- st.data:

  SpaMTP Seurat Object containing ST data

- msi.pixel.multiplier:

  Numeric value defining a scale.factor to multiple each SM pixel
  coordinates by (default = 20).

- image.res:

  Character string of the corresponding ST image scale factor to use
  (default = "lowres").

- continous_cols:

  Vector of colours to use for plotting continuous data. If NULL, the
  colour map "Reds" will be used (default = NULL).

- catagorical_cols:

  Vector of colours to use for plotting categorical data (default =
  NULL).

- fov:

  Character string matching the name of the SM FOV to use for plotting
  (default = "fov").

- image.slice:

  Character string matching the ST image slice name to use for plotting
  (default = "slice1").

- shiny.host:

  Character string of the shiny host network interface that the Shiny
  application will listen on when run (default = "0.0.0.0").

- shiny.port:

  Numeric 4 digit number defining the port that the Shiny application
  will listen to when run (default = 4698).

- verbose:

  Boolean indicating whether to show the message. If TRUE the message
  will be show, else the message will be suppressed (default = FALSE).

## Value

A SpaMTP Seurat Object containing SM data with transformed coordinated
to match the aligned ST data

## Examples

``` r
# SM_Transformed <- AlignSpatialOmics(SM.data, ST.data)
```
