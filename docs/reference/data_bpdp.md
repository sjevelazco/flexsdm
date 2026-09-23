# Calculate data to construct partial dependence surface plots

Calculate data to construct Partial dependence surface plot (i.e.,
bivariate dependence plot) for two predictor set

## Usage

``` r
data_bpdp(
  model,
  predictors,
  resolution = 50,
  training_data = NULL,
  training_boundaries = NULL,
  projection_data = NULL,
  clamping = FALSE
)
```

## Arguments

- model:

  A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm",
  "maxnet”, “nnet", and "randomForest" This model can be found in the
  first element of the list returned by any function from the fit\_,
  tune\_, or esm\_ function families

- predictors:

  character. Vector with two predictor name(s) to plot. If NULL all
  predictors will be plotted. Default NULL

- resolution:

  numeric. Number of equally spaced points at which to predict
  continuous predictors. Default 50

- training_data:

  data.frame. Database with response (0,1) and predictor values used to
  fit a model. Default NULL

- training_boundaries:

  character. Plot training conditions boundaries based on training data
  (i.e., presences, presences and absences, etc). If training_boundaries
  = "convexh", function will delimit training environmental region based
  on a convex-hull. If training_boundaries = "rectangle", function will
  delimit training environmental region based on four straight lines. If
  used any methods it is necessary provide data in training_data
  argument. If NULL all predictors will be used. Default NULL.

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. Default NULL

- clamping:

  logical. Perform clamping. Only for maxent models. Default FALSE

## Value

A list with two tibbles "pdpdata" and "resid".

- pspdata: has data to construct partial dependence surface plot, the
  first two column includes values of the selected environmental
  variables, the third column with predicted suitability.

- training_boundaries: has data to plot boundaries of training data.

## See also

[`data_pdp`](https://sjevelazco.github.io/flexsdm/reference/data_pdp.md)`, `[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)`, `[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md)

## Examples

``` r
# \donttest{
library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
data(abies)

abies2 <- abies %>%
  select(x, y, pr_ab)

abies2 <- sdm_extract(abies2,
  x = "x",
  y = "y",
  env_layer = somevar
)
#> 60 rows were excluded from database because NAs were found
abies2 <- part_random(abies2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

m <- fit_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  partition = ".part",
  thr = c("max_sens_spec")
)
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmx + tmn
#> Replica number: 1/1
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5

df <- data_bpdp(
  model = m$model,
  predictors = c("aet", "cwd"),
  resolution = 50,
  projection_data = somevar,
  training_boundaries = "rectangle",
  training_data = abies2,
  clamping = TRUE
)

df
#> $pspdata
#> # A tibble: 2,500 × 3
#>      aet   cwd Suitability
#>    <dbl> <dbl>       <dbl>
#>  1   0   -9.39      0.160 
#>  2  27.7 -9.39      0.152 
#>  3  55.4 -9.39      0.139 
#>  4  83.1 -9.39      0.121 
#>  5 111.  -9.39      0.100 
#>  6 139.  -9.39      0.0791
#>  7 166.  -9.39      0.0605
#>  8 194.  -9.39      0.0457
#>  9 222.  -9.39      0.0352
#> 10 249.  -9.39      0.0285
#> # ℹ 2,490 more rows
#> 
#> $training_boundaries
#> # A tibble: 4 × 2
#>     aet   cwd
#>   <dbl> <dbl>
#> 1  117. -6.59
#> 2 1201. -6.59
#> 3  117. 10.5 
#> 4 1201. 10.5 
#> 
names(df)
#> [1] "pspdata"             "training_boundaries"
df$pspdata
#> # A tibble: 2,500 × 3
#>      aet   cwd Suitability
#>    <dbl> <dbl>       <dbl>
#>  1   0   -9.39      0.160 
#>  2  27.7 -9.39      0.152 
#>  3  55.4 -9.39      0.139 
#>  4  83.1 -9.39      0.121 
#>  5 111.  -9.39      0.100 
#>  6 139.  -9.39      0.0791
#>  7 166.  -9.39      0.0605
#>  8 194.  -9.39      0.0457
#>  9 222.  -9.39      0.0352
#> 10 249.  -9.39      0.0285
#> # ℹ 2,490 more rows
df$training_boundaries
#> # A tibble: 4 × 2
#>     aet   cwd
#>   <dbl> <dbl>
#> 1  117. -6.59
#> 2 1201. -6.59
#> 3  117. 10.5 
#> 4 1201. 10.5 

# see p_bpdp to construct partial dependence plot with ggplot2
# }
```
