# Partial Dependent Plot

Create partial dependence plot(s) to explore the marginal effect of
predictors on suitability

## Usage

``` r
p_pdp(
  model,
  predictors = NULL,
  resolution = 100,
  resid = FALSE,
  training_data = NULL,
  projection_data = NULL,
  clamping = FALSE,
  rug = FALSE,
  colorl = c("#462777", "#6DCC57"),
  colorp = "black",
  alpha = 0.2,
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

  character. Vector of predictor name(s) to calculate partial dependence
  plots. If NULL all predictors will be used. Default NULL

- resolution:

  numeric. Number of equally spaced points at which to predict
  suitability values for continuous predictors. Default 50

- resid:

  logical. Calculate residuals based on training data. Default FALSE

- training_data:

  data.frame. Database with response (0,1) and predictor values used to
  fit a model. Default NULL

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. When this argument is used, function will calculate
  partial dependence curves distinguishing conditions used in training
  and projection conditions (i.e., projection data present in projection
  area but not training). Default NULL

- clamping:

  logical. Perform clamping. Only for maxent models. Default FALSE

- rug:

  logical. Display training data as a rug plot on the x-axis. Note: this
  could be time-consuming for large databases. Default FALSE

- colorl:

  character. A vector with one or two colors used to color lines. If
  projection_data argument is used it is necessary to provide two
  colors. Default c("#462777", "#6DCC57")

- colorp:

  character. A vector with one color used to color points of residuals,
  Default "black"

- alpha:

  numeric. a value between 0 to 1 to control transparency of residual
  points. Lower values corresponding to more transparent colors. Default
  0.2

- theme:

  ggplot2 theme. Default ggplot2::theme_classic()

## Details

This function creates partial dependent plots to explore the marginal
effect of predictors on suitability. If projection_data is used,
function will extract the minimum and maximum values found in a region
or time period to which a model will be projected. If the range of
projection data is greater than that of the training data it will be
plotted with a different color. Partial dependence curves could be used
to interpret a model or to explore how a model may extrapolate outside
the environmental conditions used to train the model.

## See also

[`data_pdp`](https://sjevelazco.github.io/flexsdm/reference/data_pdp.md),
[`data_bpdp`](https://sjevelazco.github.io/flexsdm/reference/data_bpdp.md),
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md),
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

# Partial depence plot
p_pdp(model = svm_t1$model, training_data = abies2)
p_pdp(model = svm_t1$model, training_data = abies2, predictors = c("aet", "cwd"))
p_pdp(model = svm_t1$model, training_data = abies2, resolution = 5)
p_pdp(model = svm_t1$model, training_data = abies2, resolution = 50)
p_pdp(model = svm_t1$model, training_data = abies2, resid = TRUE)
p_pdp(
  model = svm_t1$model, training_data = abies2, resid = TRUE,
  colorl = "black", colorp = "red", alpha = 0.1
)
p_pdp(
  model = svm_t1$model, training_data = abies2, resid = TRUE,
  colorl = "black", colorp = "red", alpha = 0.1, rug = TRUE
)

# Partial depence plot for training and projection condition found in a projection area
plot(somevar[[1]], main = "Projection area")
p_pdp(model = svm_t1$model, training_data = abies2, projection_data = somevar)
p_pdp(
  model = svm_t1$model, training_data = abies2, projection_data = somevar,
  colorl = c("#CC00FF", "#CCFF00")
)
p_pdp(
  model = svm_t1$model, training_data = abies2, projection_data = somevar,
  colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
)
p_pdp(
  model = svm_t1$model, training_data = abies2, projection_data = somevar,
  colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
  theme = ggplot2::theme_dark()
)
} # }
```
