# Raster interpolation between two time periods

This function interpolates values for each year between two specified
years with simple interpolation using two raster objects containing e.g.
habitat suitability values predicted using a species distribution model.

## Usage

``` r
interp(r1, r2, y1, y2, rastername = NULL, dir_save = NULL)
```

## Arguments

- r1:

  SpatRaster. Raster object for the initial year

- r2:

  SpatRaster. Raster object for the final year

- y1:

  numeric. Initial year

- y2:

  numeric. Final year

- rastername:

  character. Word used as prefix in raster file name. Default NULL

- dir_save:

  character. Directory path and name of the folder in which the raster
  files will be saved. If NULL, function will return a SpatRaster
  object, else, it will save raster in a given directory. Default NULL

## Value

If dir_save is NULL, the function returns a SpatRaster with suitability
interpolation for each year. If dir_save is used, function outputs are
saved in the directory specified in dir_save.

## Details

This function interpolates suitability values assuming that annual
changes in suitability are linear. This function could be useful for
linking SDM output based on averaged climate data and climate change
scenarios to other models that require suitability values disaggregated
in time periods, such as population dynamics (Keith et al., 2008;
Conlisk et al., 2013; Syphard et al., 2013).

## References

- Keith, D.A., Akçakaya, H.R., Thuiller, W., Midgley, G.F., Pearson,
  R.G., Phillips, S.J., Regan, H.M., Araujo, M.B. & Rebelo, T.G. (2008)
  Predicting extinction risks under climate change: coupling stochastic
  population models with dynamic bioclimatic habitat models. Biology
  Letters, 4, 560-563.

- Conlisk, E., Syphard, A.D., Franklin, J., Flint, L., Flint, A. &
  Regan, H.M. (2013) Management implications of uncertainty in assessing
  impacts of multiple landscape-scale threats to species persistence
  using a linked modeling approach. Global Change Biology 3, 858-869.

- Syphard, A.D., Regan, H.M., Franklin, J. & Swab, R. (2013) Does
  functional type vulnerability to multiple threats depend on spatial
  context in Mediterranean-climate regions? Diversity and Distributions,
  19, 1263-1274.

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)

f <- system.file("external/suit_time_step.tif", package = "flexsdm")
abma <- terra::rast(f)
plot(abma)


int <- interp(
  r1 = abma[[1]],
  r2 = abma[[2]],
  y1 = 2010,
  y2 = 2020,
  rastername = "Abies",
  dir_save = NULL
)

int
#> class       : SpatRaster
#> size        : 558, 394, 11  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> varnames    : suit_time_step
#>               suit_time_step
#>               suit_time_step
#>               suit_time_step
#>               suit_time_step
#>               ...
#> names       : Abies_2010, Abies_2011, Abies_2012, Abies_2013, Abies_2014, Abies_2015, ...
#> min values  :          0,          0,          0,          0,          0,          0, ...
#> max values  :   0.975611,   0.960608,   0.950461,   0.944007,   0.944294,   0.946355, ...
# }
```
