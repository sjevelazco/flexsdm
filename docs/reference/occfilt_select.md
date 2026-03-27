# Select filtered occurrences

Select filtered occurrences based on number of records and spatial
autocorrelation (see details)

## Usage

``` r
occfilt_select(occ_list, x, y, env_layer, filter_prop = FALSE)
```

## Arguments

- occ_list:

  list. A list with filtered specie occurrences testing several values
  (see
  [`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md)
  and
  [`occfilt_geo`](https://sjevelazco.github.io/flexsdm/reference/occfilt_geo.md))

- x:

  character. Column name with longitude data

- y:

  character. Column name with latitude data

- env_layer:

  SpatRaster. Raster variables that will be used to fit the model.
  Factor variables will be removed.

- filter_prop:

  logical. If TRUE, the function will return a list with the filtered
  occurrences and a tibble with the spatial autocorrelation and number
  of occurrence values

## Value

If filter_prop = FALSE, a tibble with selected filtered occurrences. If
filter_prop = TRUE, a list with following objects:

- A tibble with selected filtered occurrences

- A tibble with filter properties with columns:

  - filt_value: values used for filtering, the value with an asterisk
    will denote the one selected

  - n_records: number of occurrence

  - mean_autocorr: mean spatial autocorrelation.

  - the remaining columns have the spatial autocorrelation values for
    each variable.

## Details

The function implement the approach used in Velazco et al. (2020) which
consists in calculating for each filtered dataset:

- 1- the number of occurrence.

- 2- the spatial autocorrelation based on Morans'I for each variable

- 3- the mean spatial autocorrelation among variables

Then function will select those dataset with average spatial
autocorrelation lower than the mean of all dataset, and from this subset
will select the one with the highest number occurrences.

If use occfilt_select cite Velazco et al. (2020) as reference.

## References

- Velazco, S. J. E., Svenning, J-C., Ribeiro, B. R., & Laureto, L. M. O.
  (2020). On opportunities and threats to conserve the phylogenetic
  diversity of Neotropical palms. Diversity and Distributions, 27,
  512–523. https://doi.org/10.1111/ddi.13215

## See also

[`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md),
[`occfilt_geo`](https://sjevelazco.github.io/flexsdm/reference/occfilt_geo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(dplyr)

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

plot(somevar)

# Species occurrences
data("spp")
spp
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

## %######################################################%##
####                  Cellsize method                   ####
## %######################################################%##
# Using cellsize method
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("cellsize", factor = c(1, 4, 8, 12, 16, 20)),
  prj = crs(somevar)
)

filtered_occ

# Select filtered occurrences based on
# number of records and spatial autocorrelation
occ_selected <- occfilt_select(
  occ_list = filtered_occ,
  x = "x",
  y = "y",
  env_layer = somevar,
  filter_prop = FALSE
)
occ_selected

occ_selected <- occfilt_select(
  occ_list = filtered_occ,
  x = "x",
  y = "y",
  env_layer = somevar,
  filter_prop = TRUE
)
occ_selected$occ

occ_selected$filter_prop
} # }
```
