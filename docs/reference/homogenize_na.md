# Homogenize cells with NAs across all layers

Homogenize cells with NAs across all layers

## Usage

``` r
homogenize_na(x)
```

## Arguments

- x:

  A SpatRaster.

## Value

a SpatRaster

## Details

Homogenize cells with NAs across all layers in a SpatRaster so that in
the resulting SpatRaster all layers have the same cells with NAa

## Examples

``` r
if (FALSE) { # \dontrun{
#' require(terra)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

somevar2 <- homogenize_na(somevar)
par(mfrow = c(2, 1))
plot(somevar$CFP_4)
plot(somevar2$CFP_4)
par(mfrow = c(1, 1))

# In somevar2 all layers have the same cells with NAs
} # }
```
