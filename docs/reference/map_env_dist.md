# Calculate environmental distance between presences and projection data

Calculate environmental distance between presences and projection data

## Usage

``` r
map_env_dist(training_data, projection_data, metric = "domain", n_cores = 1)
```

## Arguments

- training_data:

  data.frame or tibble with environmental conditions of presence used
  for constructing models

- projection_data:

  SpatRaster, data.frame or tibble with environmental condition used for
  projecting a model (e.g., a larger, encompassing region, a spatially
  separate region, or a different time period). If data.frame or tibble
  is used function will return a tibble object. Otherwise, as SpatRaster
  object.

- metric:

  character. Metric used for measuring distance. Default 'domain'.

  - domain: Gower distance between presences and projection data. The
    result is a similarity value that ranges from 0 to 1.

  - euclidean: Euclidean distance. This metric Z-score standardizes the
    projection data based on the training data and then standardizes the
    resulting distance by the maximum pairwise distance within the
    training data. The result is a similarity value that ranges from 0
    to 1.

  - mahalanobis: Mahalanobis distance. This metric Z-score standardizes
    the projection data based on the training data and then standardizes
    the resulting distance by the maximum pairwise distance within the
    training data. The result is a similarity value that ranges from 0
    to 1.

- n_cores:

  numeric. Number of cores to use for parallel processing when metric is
  "domain". Default 1 (no parallelization).

## Value

A SpatRaster or tibble object with the nearest environmental distance
between presences and projection data. So far the Domain algorithm
(based on the Gower distance; Carpenter et al., 1993), Euclidean
distance, and Mahalanobis distance have been implemented. For Euclidean
and Mahalanobis distances, the result is a similarity value ranging from
0 to 1, where values close to 1 represent high similarity (low distance)
and values close to 0 represent low similarity (high distance).

## References

- Carpenter, G., Gillison, A.N., Winter, J., 1993. DOMAIN: a flexible
  modelling procedure for mapping potential distributions of plants and
  animals. Biodiversity & Conservation 2, 667–680

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)
data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Let's use only two variables to turn more evident the pater in the environmental space
somevar <- somevar[[1:2]]
names(somevar) <- c("aet", "cwd")


spp$species %>% unique()
sp <- spp %>%
  dplyr::filter(species == "sp3", pr_ab == 1) %>%
  dplyr::select(x, y, pr_ab)

# Get environmental condition of presences
sp_pa_2 <- sdm_extract(
  data = sp,
  x = "x",
  y = "y",
  env_layer = somevar
)
sp_pa_2

# Measure environmental distance between presences and projection data
clrs = c("#000033", "#1400FF", "#C729D6", "#FF9C63", "#FFFF60")
# Domain
envdist <-
  map_env_dist(
    training_data = sp_pa_2,
    projection_data = somevar,
    metric = "domain"
  )
plot(envdist, main = "Domain")
p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = envdist,
  projection_data = somevar,
  geo_space = FALSE,
  prop_points = 0.8,
  alpha_p = 0.9,
  color_p = "red",
  color_gradient = clrs
)
# Euclidean
envdist <-
  map_env_dist(
    training_data = sp_pa_2,
    projection_data = somevar,
    metric = "euclidean"
  )
p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = envdist,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.8,
  alpha_p = 0.9,
  color_p = "red",
  color_gradient = clrs
) +
  labs(title = "Euclidean")


# Mahalanobis
envdist <-
  map_env_dist(
    training_data = sp_pa_2,
    projection_data = somevar,
    metric = "mahalanobis"
  )
p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = envdist,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.8,
  alpha_p = 0.9,
  color_p = "red",
  color_gradient = clrs
) +
  labs(title = "Mahalanobis")
} # }
```
