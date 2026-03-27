# Plot different resolutions to be used in part_sblock

This function is useful to display the maximum and minimum resolution
that you want to test with the block_partition function. Note that if
the resolution to be tested is very fine, the plot display may take a
long time.

## Usage

``` r
plot_res(r, res_mult)
```

## Arguments

- r:

  SpatRaster. A raster layer, preferably a layer of environmental
  variables to be used

- res_mult:

  numeric. Maximum or minimum resolution to be tested.

## Value

A plot with the original raster overlapped by a grid with the resolution
used

## Examples

``` r
if (FALSE) { # \dontrun{
f <- system.file("external/somevar.tif", package = "flexsdm")
r <- terra::rast(f)
r <- r$CFP_1
plot_res(r, res_mult = 100)
plot_res(r, res_mult = 200)
} # }
```
