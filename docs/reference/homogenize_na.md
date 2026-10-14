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
somevar <- terra::rast(system.file("external/somevar.tif", package = "flexsdm"))
somevar <- terra::aggregate(somevar, fact = 6) # coarser resolution to keep this example fast

somevar2 <- homogenize_na(somevar)

# In somevar2 all layers have the same cells with NAs
terra::global(somevar2, "notNA")
#>       notNA
#> CFP_1  1488
#> CFP_2  1488
#> CFP_3  1488
#> CFP_4  1488

# \donttest{
require(terra)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

somevar2 <- homogenize_na(somevar)
oldpar <- par(mfrow = c(2, 1))
plot(somevar$CFP_4)
plot(somevar2$CFP_4)

par(oldpar) # restore graphical parameters

# In somevar2 all layers have the same cells with NAs
# }
```
