# Bivariate partial dependence plot

Create bivariate partial dependence plot(s) to explore the bivariate
marginal effect of predictors on suitability

## Usage

``` r
p_bpdp(
  model,
  predictors = NULL,
  resolution = 50,
  training_data = NULL,
  training_boundaries = NULL,
  projection_data = NULL,
  clamping = FALSE,
  color_gradient = c("#000004", "#1B0A40", "#4A0C69", "#781B6C", "#A42C5F", "#CD4345",
    "#EC6824", "#FA990B", "#F7CF3D", "#FCFFA4"),
  color_training_boundaries = "white",
  theme = ggplot2::theme_classic()
)
```

## Arguments

- model:

  A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm",
  "maxnet”, “nnet", and "randomForest" This model can be found in the
  first element of the list returned by any function from the fit\_,
  tune\_, or esm\_ function families

- predictors:

  character. Vector of predictor names to calculate partial dependence
  plots. If NULL all predictors will be used. Default NULL

- resolution:

  numeric. Number of equally spaced points at which to predict
  suitability values for continuous predictors. Default 50

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

- color_gradient:

  character. A vector with range of colors to plot. Default c("#FDE725",
  "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D",
  "#3D4988", "#462777", "#440154")

- color_training_boundaries:

  character. A vector with one color used to color points of residuals,
  Default "white"

- theme:

  ggplot2 theme. Default ggplot2::theme_classic()

## Details

This function creates partial dependent surface plots to explore the
bivariate marginal effect of predictors on suitability. If
projection_data is used, function will extract the minimum and maximum
values found in a region or time period to which a model will be
projected. Partial dependence surface plot could be used to interpret a
model or to explore how a model my extrapolate outside the environmental
conditions used to train the model (convex hull polygon).

## See also

[`data_pdp`](https://sjevelazco.github.io/flexsdm/reference/data_pdp.md),
[`data_bpdp`](https://sjevelazco.github.io/flexsdm/reference/data_bpdp.md),
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md),
[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md),
[`extra_truncate`](https://sjevelazco.github.io/flexsdm/reference/extra_truncate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
data(abies)

# set seed
abies2 <- abies %>%
  dplyr::select(x, y, pr_ab) %>%
  dplyr::group_by(pr_ab) %>%
  dplyr::slice_sample(prop = 0.5)

abies2 <- sdm_extract(abies2,
  x = "x",
  y = "y",
  env_layer = somevar
)
abies2 <- part_random(abies2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

svm_t1 <- fit_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  partition = ".part",
  thr = c("max_sens_spec")
)

# Partial depence surface plot
p_bpdp(model = svm_t1$model, training_data = abies2)
p_bpdp(model = svm_t1$model, training_data = abies2, predictors = c("aet", "cwd"))
p_bpdp(model = svm_t1$model, training_data = abies2, resolution = 10)
p_bpdp(model = svm_t1$model, training_data = abies2, resolution = 70)
# With training condition boundaires
p_bpdp(
  model = svm_t1$model, training_data = abies2,
  training_boundaries = "convexh"
)
p_bpdp(
  model = svm_t1$model, training_data = abies2,
  training_boundaries = "rectangle", color_training_boundaries = "yellow"
)
p_bpdp(
  model = svm_t1$model, training_data = abies2, training_boundaries = "convexh",
  color_training_boundaries = "orange",
  color_gradient = c("#00007F", "#007FFF", "#7FFF7F", "#FF7F00", "#7F0000")
)
# With projection data
p_bpdp(
  model = svm_t1$model, training_data = abies2, training_boundaries = "rectangle",
  projection_data = somevar, # a SpatRaster used to predict or project the model
  color_training_boundaries = "white",
  color_gradient = c("#00007F", "#007FFF", "#7FFF7F", "#FF7F00", "#7F0000")
)

# Bivariate partial dependence plot for training and projection condition
plot(somevar[[1]], main = "Projection area")
p_bpdp(
  model = svm_t1$model, training_data = abies2,
  projection_data = somevar, # a SpatRaster used to predict or project the model
  training_boundaries = "convexh"
)


# Bivariate partial dependece plot with categorical variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field = "clusters")
somevar <- c(somevar, cat)
plot(somevar)

# set seed
abies2 <- abies %>%
  dplyr::select(x, y, pr_ab) %>%
  dplyr::group_by(pr_ab) %>%
  dplyr::slice_sample(prop = 0.5)

abies2 <- sdm_extract(
  data = abies2,
  x = "x",
  y = "y",
  env_layer = somevar
)
abies2 <- part_random(abies2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

svm_t1 <- fit_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  predictors_f = "clusters",
  partition = ".part",
  thr = c("max_sens_spec")
)

p_bpdp(model = svm_t1$model, training_data = abies2, training_boundaries = "convexh")
} # }
```
