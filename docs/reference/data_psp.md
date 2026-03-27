# Calculate data to construct partial dependence surface plots

Calculate data to construct Partial dependence surface plot (i.e.,
bivariate dependence plot) for two predictor set

## Usage

``` r
data_psp(
  model,
  predictors,
  resolution = 50,
  training_data = NULL,
  pchull = FALSE,
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

- pchull:

  logical. Extract convex-hull limit of training data. Default FALSE

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

- pchull: has data to plot residuals a convex hull polygon bounding
  calibration data.

## See also

[`data_pdp`](https://sjevelazco.github.io/flexsdm/reference/data_pdp.md)`, `[`p_psp`](https://sjevelazco.github.io/flexsdm/reference/p_psp.md)`, `[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md)

## Examples

``` r
if (FALSE) {
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

df <- data_psp(
  model = m$model,
  predictors = c("aet", "cwd"),
  resolution = 50,
  projection_data = somevar,
  pchull = TRUE,
  training_data = abies2,
  clamping = TRUE
)

df
names(df)
df$pspdata
df$pchull

# see p_psp to construct partial dependence plot with ggplot2
}
```
