# Spatial band cross-validation

This function explores different numbers of spatial bands and returns
the most suitable value for a given presence or presence-absence
database. The selection of the best number of bands is performed
automatically considering spatial autocorrelation, environmental
similarity, and the number of presence and absence records in each
partition.

## Usage

``` r
part_sband(
  env_layer,
  data,
  x,
  y,
  pr_ab,
  type = "lon",
  n_part = 2,
  min_bands = 2,
  max_bands = 20,
  min_occ = 10,
  prop = 0.5
)
```

## Arguments

- env_layer:

  SpatRaster. Raster with environmental variable. Used to evaluate
  spatial autocorrelation and environmental similarity between training
  and testing partitions. Because this function calculate dissimilarity
  based on Euclidean distances, it can only be used with continuous
  environmental variables

- data:

  data.frame. Data.frame or tibble object with presences (or
  presence-absence, or presence-pseudo-absence) records, and coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- pr_ab:

  character. Column with presences, presence-absence, or
  -pseudo-absence. Presences must be represented by 1 and absences by 0

- type:

  character. Specify bands across different degrees of longitude 'lon'
  or latitude 'lat'. Default is 'lon'.

- n_part:

  integer. Number of partition. Default 2, values other than 2 has not
  yet been implemented.

- min_bands:

  integer. Minimum number of spatial bands to be tested, default 2.

- max_bands:

  integer. Maximum number of spatial bands to be tested, default 20.

- min_occ:

  numeric. Minimum number of presences or absences in a partition fold.
  The min_occ value should be base on the number of predictors in order
  to avoid over-fitting or error when fitting models for a given fold.
  Default 10.

- prop:

  numeric. Proportion of points used for testing autocorrelation between
  groups (values \> 0 and \<=1). The smaller this number is, the faster
  the function will work. Default 0.5

## Value

A list with:

- part: A tibble object with information used in 'data' arguments and a
  additional column .part with partition group.

- best_part_info: A tibble with information about the best partition. It
  contains the number of the best partition (n_grid), number of bands
  (n_bands), standard deviation of presences (sd_p), standard deviation
  of absences (sd_a), Moran's I spatial autocorrelation (spa_auto), and
  environmental similarity based on Euclidean distance (env_sim).

- grid: A SpatRaster object with bands

## Details

The part_sbands function allows testing different numbers of partitions
using a range of latitudinal or longitudinal bands. This function
explores a range of numbers of bands for a given number of partitions
and automatically selects the best number of bands for a given presence,
presence-absences, or presence-pseudo-absences dataset. Selection of
number of bands is based on an optimization procedure that explores
partitions in three dimensions determined by spatial autocorrelation
(measured by Moran's I), environmental similarity (Euclidean distance),
and difference in the amount of data among partition groups (Standard
Deviation - SD; Velazco et al., 2019). This procedure is iterative; it
will first select those partitions with autocorrelation values less than
the lowest quartile of Morans I, then those with environmental
similarity values greater than the third quartile of the Euclidean
distances, then those with a difference in the amount of data less than
the lowest quartile of SD. This selection is repeated until only one
partition is retained (Velazco et al., 2019). The main benefits of this
partition selection are that it i) is not subjective, ii) balances the
environmental similarity and special autocorrelation between partitions
groups, and iii) controls the selection of partitions with very little
data that may be problematic for model fitting ("min_occ" argument).

Partitions that are geographically structured tend to evaluate model
transferability more directly than conventional ones (e.g., those
performed by
[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md))
(Roberts et al., 2017; Santini et al., 2021), being relevant for models
that are to be used for projections in other regions outside the
training area or for other time periods. Band partitions can be an
option for those species where no best partition is found with
part_sblock or for species that are distributed linearly (e.g., species
that inhabit coastlines).

This function can interact with
[`get_block`](https://sjevelazco.github.io/flexsdm/reference/get_block.md),
[`sample_background`](https://sjevelazco.github.io/flexsdm/reference/sample_background.md),
and
[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md)
for sampling background points or pseudo-absences within spatial
partition groups

## References

- Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
  Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schroder,
  B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., &
  Dormann, C. F. (2017). Cross-validation strategies for data with
  temporal, spatial, hierarchical, or phylogenetic structure. Ecography,
  40, 913-929. https://doi.org/10.1111/ecog.02881

- Santini, L., Benitez-Lopez, A., Maiorano, L., Cengic, M., &
  Huijbregts, M. A. J. (2021). Assessing the reliability of species
  distribution projections in climate change research. Diversity and
  Distributions, ddi.13252. https://doi.org/10.1111/ddi.13252

- Velazco, S. J. E., Villalobos, F., Galvao, F., & De Marco Junior, P.
  (2019). A dark scenario for Cerrado plant species: Effects of future
  climate, land use and protected areas ineffectiveness. Diversity and
  Distributions, 25(4), 660-673. https://doi.org/10.1111/ddi.12886

## See also

[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md),
[`part_sblock`](https://sjevelazco.github.io/flexsdm/reference/part_sblock.md),
[`part_senv`](https://sjevelazco.github.io/flexsdm/reference/part_senv.md),
and
[`get_block`](https://sjevelazco.github.io/flexsdm/reference/get_block.md)

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)

# Load datasets
data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Example of two longitudinal partitions with presences and absences
single_spp <- spp %>% dplyr::filter(species == "sp1")
part_1 <- part_sband(
  env_layer = somevar,
  data = single_spp,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  type = "lon",
  min_bands = 2,
  max_bands = 20,
  n_part = 2,
  min_occ = 10,
  prop = 0.5
)
#> 55 rows were excluded from database because NAs were found
#> The following number of bands will be tested:
#> 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20
#> Creating basic raster mask...
#> Searching for the optimal number of bands...

part_1$part # database with partition fold (.part)
#> # A tibble: 945 × 4
#>           x        y pr_ab .part
#>       <dbl>    <dbl> <dbl> <int>
#>  1   -5541. -145138.     0     2
#>  2  -51981.   16322.     0     2
#>  3 -269871.   69512.     1     1
#>  4  -96261.  -32008.     0     2
#>  5  269589. -566338.     0     1
#>  6   29829. -328468.     0     2
#>  7 -152691.  393782.     0     1
#>  8 -195081.  253652.     0     1
#>  9    -951. -277978.     0     2
#> 10  145929. -271498.     0     1
#> # ℹ 935 more rows
part_1$part %>%
  group_by(pr_ab, .part) %>%
  count() # number of presences and absences in each fold
#> # A tibble: 4 × 3
#> # Groups:   pr_ab, .part [4]
#>   pr_ab .part     n
#>   <dbl> <int> <int>
#> 1     0     1   345
#> 2     0     2   362
#> 3     1     1   171
#> 4     1     2    67
part_1$best_part_info # information of the best partition
#> # A tibble: 1 × 6
#>   n_grid n_bands spa_auto env_sim  sd_p  sd_a
#>    <int>   <int>    <dbl>   <dbl> <dbl> <dbl>
#> 1      2       3    0.191    356.  73.5  12.0
part_1$grid # raster with folds
#> class       : SpatRaster
#> size        : 1, 3, 1  (nrow, ncol, nlyr)
#> resolution  : 248220, 1054620  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> name        : .part
#> min value   :     1
#> max value   :     2

# Explore grid object and presences and absences points
plot(part_1$grid, col = gray.colors(20))
points(part_1$part[c("x", "y")],
  col = rainbow(8)[part_1$part$.part],
  cex = 0.9,
  pch = c(1, 19)[part_1$part$pr_ab + 1]
)



# Example of four latitudinal partition and only presences
single_spp <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
part_2 <- part_sband(
  env_layer = somevar,
  data = single_spp,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  type = "lat",
  min_bands = 8,
  max_bands = 40,
  n_part = 8,
  min_occ = 10,
  prop = 0.5
)
#> 12 rows were excluded from database because NAs were found
#> The following number of bands will be tested:
#> 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40
#> Creating basic raster mask...
#> Searching for the optimal number of bands...

part_2$part
#> # A tibble: 238 × 4
#>           x        y pr_ab .part
#>       <dbl>    <dbl> <dbl> <int>
#>  1 -269871.   69512.     1     7
#>  2 -149991.  267962.     1     4
#>  3 -126231.  196142.     1     5
#>  4   91659. -156748.     1     3
#>  5 -210471.  326282.     1     3
#>  6 -140541.  284972.     1     3
#>  7 -217491.   65732.     1     7
#>  8 -201831.   17132.     1     8
#>  9  -40101.  150782.     1     6
#> 10 -178611.  225032.     1     4
#> # ℹ 228 more rows
part_2$best_part_info
#> # A tibble: 1 × 5
#>   n_grid n_bands spa_auto env_sim  sd_p
#>    <int>   <int>    <dbl>   <dbl> <dbl>
#> 1      6      18    0.632    272.  10.1
part_2$grid
#> class       : SpatRaster
#> size        : 18, 1, 1  (nrow, ncol, nlyr)
#> resolution  : 744660, 58590  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> name        : .part
#> min value   :     1
#> max value   :     8

# Explore Grid object and presences points
plot(part_2$grid, col = gray.colors(20))
points(part_2$part[c("x", "y")],
  col = rainbow(8)[part_2$part$.part],
  cex = 0.5,
  pch = 19
)

# }
```
