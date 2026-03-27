# Transform a spatial partition layer to the same spatial properties as environmental variables

Transform a spatial partition layer to the same spatial properties as
environmental variables

## Usage

``` r
get_block(env_layer, best_grid)
```

## Arguments

- env_layer:

  SpatRaster object with some environmental variables used in the
  block_partition or band_partition function. Function always will
  select the first layer

- best_grid:

  SpatRaster object returned by block_partition or band_partition

## Value

A SpatRaster layer with the same resolution and extent as the
environmental variables

## Details

Transform a layer originating from the function block_partition or
band_partition to the same spatial properties as the environmental
variables

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)
data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Example for a single species
single_spp <- spp %>% dplyr::filter(species == "sp3")

part <- part_sblock(
  env_layer = somevar,
  data = single_spp,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  min_res_mult = 100,
  max_res_mult = 500,
  num_grids = 10,
  min_occ = 5,
  n_part = 2
)

grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
grid_env
part$grid

plot(part$grid)
plot(grid_env)
} # }
```
