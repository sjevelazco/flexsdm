# Environmental and spatial cross-validation

This function explores different numbers of environmental partitions
(clusters) based on the K-means clustering algorithm and returns the
number of partitions best suited for a given presence,
presence-absences, or presence-pseudo-absences database. Selection of
the best number of partitions is performed automatically considering
spatial autocorrelation, environmental similarity, and the number of
presence and/or absence records in each partition.

## Usage

``` r
part_senv(
  env_layer,
  data,
  x,
  y,
  pr_ab,
  min_n_groups = 2,
  max_n_groups = 10,
  min_occ = 10,
  prop = 0.5,
  include_coords = TRUE
)
```

## Arguments

- env_layer:

  SpatRaster. Raster with environmental variable. This will be used to
  evaluate spatial autocorrelation and environmental similarity between
  training and testing partitions. Because this function calculate
  dissimilarity based on Euclidean distances, it can only be used with
  continuous variables

- data:

  data.frame. Data.frame or tibble object with presence (or
  presence-absence, or presences-pseudo-absence) records, and
  coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- pr_ab:

  character. Column with presences, presence-absence, or pseudo-absence.
  Presences must be represented by 1 and absences by 0

- min_n_groups:

  integer. Minimum number of groups to be tested. Default 2.

- max_n_groups:

  integer. Maximum number of groups to be tested. Default 10.

- min_occ:

  numeric. Minimum number of presences or absences in a partition fold.
  The min_occ value should be base on the amount of predictors in order
  to avoid over-fitting or error when fitting models for a given fold.
  Default 10.

- prop:

  numeric. Proportion of point used for testing autocorrelation between
  groups (values \> 0 and \<=1). The smaller this number is, the faster
  the function will work. Default 0.5

- include_coords:

  logical. Whether to include coordinates in the clustering. If False,
  only environmental variables will be used for clustering. Default TRUE

## Value

A list with:

- part: A tibble object with information used in 'data' arguments and a
  additional column .part with partition group.

- best_part_info: A tibble with information about the best partition. It
  contains the number of partition (n_groups), standard deviation of
  presences (sd_p), standard deviation of absences (sd_a), Moran's I
  spatial autocorrelation (spa_auto) and environmental similarity based
  on Euclidean distance (env_sim)

## Details

The part_sblock allows test with different numbers of partitions defined
in the envirnomental clusters delimited the K-mean cluster algorithm.
This function explores a range of environmental clusters and
automatically selects best number of cluster for a given given presence,
presence-absences, or presence-pseudo-absences dataset. Such selection
of number of clusters is based on an optimization procedure that
explores partition size in three dimensions determined by spatial
autocorrelation (measured by Moran's I), environmental similarity
(Euclidean distance), and difference in the amount of data among
clusters (Standard Deviation - SD; Velazco et al., 2019). This procedure
will cyclically select those partitions with autocorrelation values less
than the lowest quartile of Morans I, then those with environmental
similarity values greater than the third quartile of the Euclidean
distances than those with a difference in the amount of data less than
the lowest quartile of SD. This selection is repeated until only one
partition is retained (Velazco et al., 2019). The main benefit of this
partition selection are i) this is not subjective, ii) balances the
environmental similarity and special autocorrelation between partitions,
and iii) controls the partition selection with few data that may be
problematic for model fitting ("min_occ" argument)..

Partitions geographically structured tend to evaluate model
transferability more directly than conventional ones (e.g., those
performed by
[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md))
(Roberts et al., 2017; Santini et al., 2021), being relevant for models
that want to be used for projections in other regions outside the
calibration area or for other periods.

This function can interact with
[`get_block`](https://sjevelazco.github.io/flexsdm/reference/get_block.md),
[`sample_background`](https://sjevelazco.github.io/flexsdm/reference/sample_background.md),
and
[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md)
for sampling background points or pseudo-absences within spatial
partition broups

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
and
[`part_sband`](https://sjevelazco.github.io/flexsdm/reference/part_sband.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(ggplot2)

f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Select a species
spp1 <- spp %>% dplyr::filter(species == "sp1")

part1 <- part_senv(
  env_layer = somevar,
  data = spp1,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  min_n_groups = 2,
  max_n_groups = 10,
  min_occ = 10,
  prop = 0.2
)

part1

ggplot(part1$part, aes(x, y, col = factor(.part))) +
  geom_point(aes(shape = factor(pr_ab)))

ggplot(part1$part, aes(x, y, col = factor(.part))) +
  geom_point(aes(shape = factor(pr_ab))) +
  facet_wrap(. ~ .part)

ggplot(part1$part, aes(x, y, col = factor(.part))) +
  geom_point(aes(shape = factor(pr_ab))) +
  facet_wrap(. ~ pr_ab)
} # }
```
