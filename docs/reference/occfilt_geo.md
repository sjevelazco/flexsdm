# Perform geographical filtering on species occurrences

This function perform geographical filtering of species occurrences
based on different approach to define the minimum nearest-neighbor
distance between points.

## Usage

``` r
occfilt_geo(
  data,
  x,
  y,
  env_layer,
  method,
  prj = "+proj=longlat +datum=WGS84",
  reps = 20
)
```

## Arguments

- data:

  data.frame. Data.frame or tibble object with presences (or
  presence-absence) records, and coordinates

- x:

  character. Column name with longitude data

- y:

  character. Column name with latitude data

- env_layer:

  SpatRaster. Raster variables that will be used to fit the model.
  Factor variables will be removed.

- method:

  character. Method to perform geographical thinning. Pairs of points
  are filtered based on a geographical distance criteria. For the three
  method, it is possible to use several values. If several values are
  provided, the function will return a list with the results. The
  following methods are available:

  - moran: records are filtered based on the smallest distance to the
    Moran's I value provided. If no Moran's I values is provided it will
    use 0.1. Usage method: method = c('moran') or method = c('moran',
    val = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35)).

  - cellsize: records are filtered based on the resolution of the
    environmental variables which can be aggregated to coarser
    resolution defined by the factor. Usage method: method =
    c('cellsize', factor = '2') or method = c('cellsize', factor = c(1,
    4, 8)).

  - defined: records are filtered based on a distance value (d) provided
    in km. Usage method: method = c('defined', d = c(20, 40, 60, 80)).

- prj:

  character. Projection string (PROJ4) for occurrences. Not necessary if
  the projection used is WGS84 ("+proj=longlat +datum=WGS84"). Default
  "+proj=longlat +datum=WGS84"

- reps:

  integer. Number of times to repeat the thinning process. Default 20

## Value

If one value is used to filter occurrence function will return a tibble
object with filtered data. If several values are used to filter
occurrences, the function will return a list of tibbles with filtered
data.

## Details

In this function three alternatives are implemented to determine the
distance threshold between pair of points: '

- "moran" determines the minimum nearest-neighbor distance that
  approximate to the spatial autocorrelation in occurrence data,
  following a Moran's I. To do so, a Principal Component Analysis with
  the environmental variables is performed and then the first Principal
  Component is used to calculate the semivariograms. Because of this,
  this method only allow the use of continuous variables. Sometimes,
  this method can (too) greatly reduce the number of presences.

- "cellsize" filters occurrences based on the predictors' resolution.
  This method will calculate the distance between the first two cells of
  the environmental variable and use this distance as minimum
  nearest-neighbor distance to filter occurrences. The resolution of the
  raster is aggregated based on the values used in "factor". Thus, the
  distance used for filtering can be adjusted to represent a larger grid
  size.

- "determined" this method uses any minimum nearest-neighbor distance
  specified in km.

The "thin" function from spThin package is used to filter data
(Aiello-Lammens et al., 2015) with the following argument settings reps
= 20, write.files = FALSE, locs.thinned.list.return = TRUE, and
write.log.file = FALSE.

In the case that more than one value was used in some of the methods, it
is possible to use use
[`occfilt_select`](https://sjevelazco.github.io/flexsdm/reference/occfilt_select.md)
function to automatically select a filtered database based on number of
records and spatial autocorrelation.

## References

- Aiello-Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., &
  Anderson, R. P. (2015). spThin: An R package for spatial thinning of
  species occurrence records for use in ecological niche models.
  Ecography, 38(5), 541-545. https://doi.org/10.1111/ecog.01132

## See also

[`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md),
[`occfilt_select`](https://sjevelazco.github.io/flexsdm/reference/occfilt_select.md)

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)
require(ggplot2)

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

plot(somevar)


# Species occurrences
data("spp")
spp
#> # A tibble: 1,150 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1       -5541. -145138.     0
#>  2 sp1      -51981.   16322.     0
#>  3 sp1     -269871.   69512.     1
#>  4 sp1      -96261.  -32008.     0
#>  5 sp1      269589. -566338.     0
#>  6 sp1       29829. -328468.     0
#>  7 sp1     -152691.  393782.     0
#>  8 sp1     -195081.  253652.     0
#>  9 sp1        -951. -277978.     0
#> 10 sp1      145929. -271498.     0
#> # ℹ 1,140 more rows
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

somevar[[1]] %>% plot()
points(spp1 %>% select(x, y))


## %######################################################%##
####                  Cellsize method                   ####
## %######################################################%##

# Using cellsize method
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("cellsize", factor = "3"),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Factor: x3
#> Distance threshold (km): 4.617
#> Number of filtered records: 212

somevar[[1]] %>% plot(col = gray.colors(10))
points(spp1 %>% select(x, y)) # raw data
points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data



# Using cellsize method with several values
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("cellsize", factor = c(1, 4, 8, 12, 16, 20)),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Factor: x1
#> Distance threshold (km): 1.539
#> Number of filtered records: 233
#> Factor: x4
#> Distance threshold (km): 6.157
#> Number of filtered records: 203
#> Factor: x8
#> Distance threshold (km): 12.313
#> Number of filtered records: 156
#> Factor: x12
#> Distance threshold (km): 18.47
#> Number of filtered records: 118
#> Factor: x16
#> Distance threshold (km): 24.626
#> Number of filtered records: 96
#> Factor: x20
#> Distance threshold (km): 30.783
#> Number of filtered records: 78

filtered_occ # Note that several values are provided for any filtering method
#> $`1`
#> # A tibble: 233 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -126231.  196142.     1
#>  4 sp1       91659. -156748.     1
#>  5 sp1     -210471.  326282.     1
#>  6 sp1     -140541.  284972.     1
#>  7 sp1     -217491.   65732.     1
#>  8 sp1     -178611.  225032.     1
#>  9 sp1      -92481.  155642.     1
#> 10 sp1     -367611.  266072.     1
#> # ℹ 223 more rows
#> 
#> $`4`
#> # A tibble: 203 × 4
#>    species        x       y pr_ab
#>    <chr>      <dbl>   <dbl> <dbl>
#>  1 sp1     -269871.  69512.     1
#>  2 sp1     -149991. 267962.     1
#>  3 sp1     -126231. 196142.     1
#>  4 sp1     -210471. 326282.     1
#>  5 sp1     -140541. 284972.     1
#>  6 sp1     -217491.  65732.     1
#>  7 sp1     -178611. 225032.     1
#>  8 sp1      -92481. 155642.     1
#>  9 sp1     -367611. 266072.     1
#> 10 sp1     -109491.  46292.     1
#> # ℹ 193 more rows
#> 
#> $`8`
#> # A tibble: 156 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1       91659. -156748.     1
#>  4 sp1     -210471.  326282.     1
#>  5 sp1     -140541.  284972.     1
#>  6 sp1     -217491.   65732.     1
#>  7 sp1      -92481.  155642.     1
#>  8 sp1     -367611.  266072.     1
#>  9 sp1     -184551.  425372.     1
#> 10 sp1     -260151.   57632.     1
#> # ℹ 146 more rows
#> 
#> $`12`
#> # A tibble: 118 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -140541.  284972.     1
#>  4 sp1     -178611.  225032.     1
#>  5 sp1      -34431.  212072.     1
#>  6 sp1     -331701.  304412.     1
#>  7 sp1     -215871.  -51178.     1
#>  8 sp1      -60351. -307138.     1
#>  9 sp1     -331161.  387842.     1
#> 10 sp1      -88971.  137282.     1
#> # ℹ 108 more rows
#> 
#> $`16`
#> # A tibble: 96 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -210471.  326282.     1
#>  2 sp1     -313611.  400802.     1
#>  3 sp1     -202911.  -64138.     1
#>  4 sp1      -34431.  212072.     1
#>  5 sp1     -172941. -113008.     1
#>  6 sp1      -60351. -307138.     1
#>  7 sp1     -304971.  219632.     1
#>  8 sp1     -195351.   85982.     1
#>  9 sp1      -23091.   56012.     1
#> 10 sp1     -114621.   54122.     1
#> # ℹ 86 more rows
#> 
#> $`20`
#> # A tibble: 78 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -126231.  196142.     1
#>  2 sp1     -367611.  266072.     1
#>  3 sp1     -109491.   46292.     1
#>  4 sp1      -34431.  212072.     1
#>  5 sp1     -143781.  233942.     1
#>  6 sp1     -215871.  -51178.     1
#>  7 sp1      -60351. -307138.     1
#>  8 sp1     -306321.  134852.     1
#>  9 sp1     -195351.   85982.     1
#> 10 sp1     -328731.  405392.     1
#> # ℹ 68 more rows
#> 
# fuction will return a list of tibbles with the results.
# So user must select the desired filtered dataset

# Let's explore the results
bind_rows(filtered_occ, .id = "cellSize") %>%
  dplyr::mutate(cellSize = as.numeric(cellSize)) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  facet_wrap(~cellSize)



## %######################################################%##
####                   Defined method                   ####
## %######################################################%##
# Using defined method
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("defined", d = "30"),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Distance threshold (km): 30
#> Number of filtered records: 78

somevar[[1]] %>% plot(col = gray.colors(10))
points(spp1 %>% select(x, y)) # raw data
points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data


# Using defined method with several values
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("defined", factor = c(5, 10, 15, 30, 35, 40)),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Distance threshold (km): 5
#> Number of filtered records: 210
#> Distance threshold (km): 10
#> Number of filtered records: 176
#> Distance threshold (km): 15
#> Number of filtered records: 141
#> Distance threshold (km): 30
#> Number of filtered records: 78
#> Distance threshold (km): 35
#> Number of filtered records: 70
#> Distance threshold (km): 40
#> Number of filtered records: 62

bind_rows(filtered_occ, .id = "cellSize") %>%
  dplyr::mutate(cellSize = as.numeric(cellSize)) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  facet_wrap(~cellSize)



## %######################################################%##
####                  Moran's I method                  ####
## %######################################################%##

# Using Moran's I method
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("moran"),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Moran's I threshold closest to the supplied value: 0.099
#> Distance threshold (km) : 334.908
#> Number of filtered records: 4

somevar[[1]] %>% plot(col = gray.colors(10))
points(spp1 %>% select(x, y)) # raw data
points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data


# Using Moran's I method with several values
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("moran", c(0.05, 0.15, 0.2, 0.5, 0.7)),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Moran's I threshold closest to the supplied value: 0.05
#> Distance threshold (km) : 362.625
#> Number of filtered records: 4
#> Moran's I threshold closest to the supplied value: 0.15
#> Distance threshold (km) : 304.881
#> Number of filtered records: 5
#> Moran's I threshold closest to the supplied value: 0.197
#> Distance threshold (km) : 279.473
#> Number of filtered records: 5
#> Moran's I threshold closest to the supplied value: 0.499
#> Distance threshold (km) : 140.887
#> Number of filtered records: 13
#> Moran's I threshold closest to the supplied value: 0.7
#> Distance threshold (km) : 60.046
#> Number of filtered records: 34

bind_rows(filtered_occ, .id = "moran") %>%
  dplyr::mutate(moran = as.numeric(moran)) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  facet_wrap(~moran)


# It is possible select the best of filtered
# datasets using the occfilt_selec function

occ_selected <- occfilt_select(
  occ_list = filtered_occ,
  x = "x",
  y = "y",
  env_layer = somevar,
  filter_prop = TRUE
)
#> Dataset with filtered value 0.7 was selected

occ_selected
#> $occ
#> # A tibble: 34 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -140541.  284972.     1
#>  2 sp1     -178611.  225032.     1
#>  3 sp1     -109491.   46292.     1
#>  4 sp1      -34431.  212072.     1
#>  5 sp1     -215871.  -51178.     1
#>  6 sp1     -172941. -113008.     1
#>  7 sp1     -154311.   95972.     1
#>  8 sp1      -80601.  143492.     1
#>  9 sp1     -318201.  314402.     1
#> 10 sp1      -19581.   98942.     1
#> # ℹ 24 more rows
#> 
#> $filter_prop
#>   filt_value mean_autocorr n_records      CFP_1      CFP_2      CFP_3     CFP_4
#> 1       0.05    0.18044575         4 0.10512829 0.29585163 0.14455848 0.1762446
#> 2       0.15    0.13171875         5 0.02985195 0.17312526 0.18938125 0.1345165
#> 3        0.2    0.18600777         5 0.16528010 0.27779062 0.08027048 0.2206899
#> 4        0.5    0.09863698        13 0.18494697 0.04271891 0.01368967 0.1531924
#> 5      * 0.7    0.12790516        34 0.17011719 0.01647243 0.12130328 0.2037277
#> 
# }
```
