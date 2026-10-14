# Calculate species distribution model uncertainty

This function calculates the uncertainty of a species distribution model
by performing a bootstrap procedure. It refits the model multiple times
on resampled data and then calculates the standard deviation of the
predictions across all iterations.

## Usage

``` r
sdm_uncertainty(
  models,
  training_data,
  background = NULL,
  response,
  projection_data,
  iteration = 50,
  n_cores = 5,
  clamp = TRUE,
  pred_type = "cloglog"
)
```

## Arguments

- models:

  A flexsdm model object from \`fit\_\*\` and \`tune\_\*\` functions.

- training_data:

  A data.frame or tibble with presence-absence/background data and
  predictors.

- background:

  A data.frame or tibble with background data, used for Maxent and
  Gausian Process models. Default NULL.

- response:

  character. Column name of the response variable.

- projection_data:

  A SpatRaster object with the environmental layers for projection.

- iteration:

  numeric. The number of bootstrap iterations. A large number of
  iterations will increase computation time, especially for large areas
  or high-resolution data. Default 50.

- n_cores:

  numeric. The number of cores to use for parallel processing. Default
  5.

- clamp:

  logical. If TRUE, predictors and features are restricted to the range
  seen during model training. Only for Maxent models. Default TRUE.

- pred_type:

  character. The type of response required. Available options are:
  "link", "exponential", "cloglog", and "logistic". Only for Maxent
  models. Default "cloglog".

## Value

A SpatRaster object with a single layer representing the model
uncertainty, calculated as the standard deviation of the bootstrap
predictions.

## See also

[`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)

## Examples

``` r
# Toy data: small raster and records of one species to keep this example fast
data("spp")
somevar <- terra::rast(system.file("external/somevar.tif", package = "flexsdm"))
somevar <- terra::aggregate(somevar, fact = 6) # coarser resolution
sp3 <- spp[spp$species == "sp3", ]
sp3 <- sdm_extract(data = sp3, x = "x", y = "y", env_layer = somevar)
#> 20 rows were excluded from database because NAs were found
sp3 <- part_random(
  data = sp3,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

glm_toy <- fit_glm(
  data = sp3,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  thr = "max_sens_spec"
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4 + I(CFP_1^2) + I(CFP_2^2) + I(CFP_3^2) + I(CFP_4^2)
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

# Few bootstrap iterations to keep this example fast
unc <- sdm_uncertainty(
  models = glm_toy,
  training_data = sp3,
  response = "pr_ab",
  projection_data = somevar,
  iteration = 2,
  n_cores = 1
)
#> 
#> Attaching package: ‘kernlab’
#> The following object is masked from ‘package:ggplot2’:
#> 
#>     alpha
#> The following objects are masked from ‘package:terra’:
#> 
#>     buffer, size
terra::plot(unc)


# \donttest{
require(terra)
require(dplyr)

data(spp)
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

sp_data <- spp %>% filter(species == "sp3")
sp_data <- part_random(
  data = sp_data,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

sp_data <- sdm_extract(
  data = sp_data,
  x = "x",
  y = "y",
  env_layer = somevar
)
#> 4 rows were excluded from database because NAs were found

m <- fit_svm(
  data = sp_data,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part"
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

unc <- sdm_uncertainty(
  models = m,
  training_data = sp_data,
  response = "pr_ab",
  projection_data = somevar,
  iteration = 10,
  n_cores = 2
)

plot(unc)

# }
```
