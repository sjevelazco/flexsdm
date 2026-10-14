# Truncate suitability predictions based on an extrapolation value

Exclusion of suitability predictions in environmental conditions assumed
to be extrapolative.

## Usage

``` r
extra_truncate(suit, extra, threshold = 50, trunc_value = 0)
```

## Arguments

- suit:

  SpatRaster with suitability values

- extra:

  SpatRaster with extrapolation values preferable measured with
  extra_eval function

- threshold:

  numeric. Vector with one or more extrapolation values used for
  truncate suitability Default 50%

- trunc_value:

  numeric. Numeric value to be used to those cells assumed to be
  extrapolative

## Value

A SpatRaster object with truncated suitability values

## Details

Exclusion of suitability predictions in environmental conditions assumed
to be extrapolative. In this function it is possible to use any metric
measuring degree of extrapolation (e.g., MESS-Multivariate Environmental
Similarity Surfaces, EO-Environmental Overlap, MOP-Mobility-Oriented
Parity, EXDET-Extrapolation Detection, or AOA-Area of Applicability).
However, we recommend to use Shape approach (see
[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md),
and Velazco et al., 2023;
[doi:10.1111/ecog.06992](https://doi.org/10.1111/ecog.06992) ).

This function truncates suitability predictions assigning a given value,
generally 0 or NA. Usage trunc_value = NA. Default 0.

To those cells assumed to be extrapolative, i.e., higher than a given
threshold of a given extrapolation metric.

See this [vignette at flexsdm
website](https://sjevelazco.github.io/flexsdm/articles/v06_Extrapolation_example.html)
for further details about Shape metric, model truncation, and tools to
explore model extrapolation.

## See also

[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md),
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md),
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md),
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)

## Examples

``` r
# Toy data: small raster and records of one species to keep this example fast
data("spp")
somevar <- terra::rast(system.file("external/somevar.tif", package = "flexsdm"))
somevar <- terra::aggregate(somevar, fact = 6) # coarser resolution
sp3 <- spp[spp$species == "sp3", ]
sp3 <- sdm_extract(data = sp3, x = "x", y = "y", env_layer = somevar)
#> 20 rows were excluded from database because NAs were found

# Extrapolation (Shape metric) in the projection area
extr <- extra_eval(
  training_data = sp3,
  pr_ab = "pr_ab",
  projection_data = somevar,
  metric = "mahalanobis"
)

# Any suitability layer can be truncated; here, a toy layer rescaled to 0-1
suit <- somevar[[1]]
suit <- suit / terra::global(suit, "max", na.rm = TRUE)[[1]]
names(suit) <- "suit"

# Suitability is set to 0 where extrapolation is higher than the threshold
suit_trunc <- extra_truncate(
  suit = suit,
  extra = extr,
  threshold = 50,
  trunc_value = 0
)
suit_trunc # one raster for each threshold value
#> $`50`
#> class       : SpatRaster
#> size        : 93, 66, 1  (nrow, ncol, nlyr)
#> resolution  : 11340, 11340  (x, y)
#> extent      : -373685.8, 374754.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> name        : suit
#> min value   :    0
#> max value   :    1
#> 

# see more examples in extra_eval function
```
