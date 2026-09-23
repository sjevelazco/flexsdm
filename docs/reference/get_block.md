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
# \donttest{
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
#> 4 rows were excluded from database because NAs were found
#> The following grid cell sizes will be tested:
#> 189000 | 273000 | 357000 | 441000 | 525000 | 609000 | 693000 | 777000 | 861000 | 945000
#> Creating basic raster mask...
#> Searching for the optimal grid size...

grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
grid_env
#> class       : SpatRaster
#> size        : 558, 394, 1  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> varname     : somevar
#> name        : .part
#> min value   :     1
#> max value   :     2
part$grid
#> class       : SpatRaster
#> size        : 6, 5, 1  (nrow, ncol, nlyr)
#> resolution  : 273000, 273000  (x, y)
#> extent      : -646685.8, 718314.2, -877813.3, 760186.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> name        : .part
#> min value   :     1
#> max value   :     2

plot(part$grid)

plot(grid_env)

# }
```
