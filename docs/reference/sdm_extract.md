# Extract environmental data values from a spatial raster based on x and y coordinates

Extract environmental data values from a spatial raster based on x and y
coordinates

## Usage

``` r
sdm_extract(data, x, y, env_layer, variables = NULL, filter_na = TRUE)
```

## Arguments

- data:

  data.frame. Database with species presence, presence-absence, or
  pseudo-absence records with x and y coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- env_layer:

  SpatRaster. Raster or raster stack with environmental variables.

- variables:

  character. Vector with the variable names of predictor (environmental)
  variables Usage variables. = c("aet", "cwd", "tmin"). If no variable
  is specified, function will return data for all layers. Default NULL

- filter_na:

  logical. If filter_na = TRUE (default), the rows with NA values for
  any of the environmental variables are removed from the returned
  tibble.

## Value

A tibble that returns the original data base with additional columns for
the extracted environmental variables at each xy location from the
SpatRaster object used in 'env_layer'

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)

# Load datasets
data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Extract environmental data from somevar for all locations in spp
ex_spp <-
  sdm_extract(
    data = spp,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = NULL,
    filter_na = FALSE
  )

# Extract environmental for two variables and remove rows with NAs
ex_spp2 <-
  sdm_extract(
    data = spp,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = c("CFP_3", "CFP_4"),
    filter_na = TRUE
  )

ex_spp
ex_spp2
} # }
```
