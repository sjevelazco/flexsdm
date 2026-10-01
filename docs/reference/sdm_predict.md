# Spatial predictions from individual and ensemble models

This function allows the geographical prediction of one or more models
constructed with the fit\_ or tune\_ function set, models fitted with
esm\_ function set (i.e., ensemble of small models approach), or models
constructed with fit_ensemble function. It can return continuous or
continuous and binary predictions for one or more thresholds

## Usage

``` r
sdm_predict(
  models,
  pred,
  nchunk = 1,
  thr = NULL,
  con_thr = FALSE,
  predict_area = NULL,
  clamp = TRUE,
  pred_type = "cloglog"
)
```

## Arguments

- models:

  list of one or more models fitted with fit\_ or tune\_ functions. In
  case use models fitted with fit_ensemble or esm\_ family function only
  one model could be used. Usage models = mglm or models = list(mglm,
  mraf, mgbm)

- pred:

  SpatRaster. Raster layer with predictor variables. Names of layers
  must exactly match those used in model fitting.

- nchunk:

  integer. Number of chunks to split data used to predict models (i.e.,
  SpatRaster used in pred argument). Predicting models in chunks helps
  reduce memory requirements in cases where models are predicted for
  large scales and high resolution. Default = 1

- thr:

  character. Threshold used to get binary suitability values (i.e.,
  0,1). It is possible to use more than one threshold type. It is
  mandatory to use the same threshold/s used to fit the models. The
  following threshold types are available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the sensitivity and specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which FPB is highest.

  - sensitivity: Threshold based on a specified sensitivity value used
    to fit the models.

  - all: All the threshold used in the model outputs used in 'models'
    argument will be used.

  Usage thr = c('lpt', 'max_sens_spec', 'max_jaccard'), thr=c('lpt',
  'max_sens_spec', 'sensitivity'), or thr='all'. If no threshold is
  specified (i.e., thr = NULL) function will return continuous
  prediction only. Default NULL

- con_thr:

  logical. If true predictions with suitability values above threshold/s
  will be returned. Default = FALSE

- predict_area:

  SpatVector, SpatialPolygon, or SpatialPolygonDataFrame. Spatial
  polygon used for restring prediction into only a given region. Default
  = NULL

- clamp:

  logical. It is set with TRUE, predictors and features are restricted
  to the range seen during model training. Only valid for Maxent model
  (see tune_mx and fit_mx). Default TRUE.

- pred_type:

  character. Type of response required available "link", "exponential",
  "cloglog" and "logistic". Only valid for Maxent model (see tune_mx and
  fit_mx). Default "cloglog".

## Value

A list of SpatRaster with continuous and/or binary predictions

## See also

[`sdm_uncertainty`](https://sjevelazco.github.io/flexsdm/reference/sdm_uncertainty.md)

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)

data("spp")
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Extract data
some_sp <- spp %>%
  filter(species == "sp3")

some_sp <-
  sdm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = somevar
  )
#> 4 rows were excluded from database because NAs were found

# Partition
some_sp <- part_random(
  data = some_sp,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 3, replicates = 5)
)


## %######################################################%##
#                                                          #
####          Create different type of models           ####
#                                                          #
## %######################################################%##
# Fit some models
mglm <- fit_glm(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  poly = 2
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4 + I(CFP_1^2) + I(CFP_2^2) + I(CFP_3^2) + I(CFP_4^2)
#> Replica number: 1/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 2/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 3/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 4/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 5/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
mraf <- fit_raf(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4
#> Replica number: 1/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 2/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 3/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 4/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 5/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
mgbm <- fit_gbm(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part"
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4
#> Replica number: 1/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 2/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 3/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 4/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Replica number: 5/5
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

# Fit an ensemble model
mensemble <- fit_ensemble(
  models = list(mglm, mraf, mgbm),
  ens_method = "meansup",
  thr = NULL,
  thr_model = "max_sens_spec",
  metric = "TSS"
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

# Fit a model with the Ensembles of Small Models approach
# Without threshold specification and with kfold
msmall <- esm_gam(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  thr = NULL
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%


## %######################################################%##
#                                                           #
####      Predict different kind of models               ####
#                                                           #
## %######################################################%##

# sdm_predict can be used for predict one or more models fitted with fit_ or tune_ functions

# a single model
ind_p <- sdm_predict(
  models = mglm,
  pred = somevar,
  thr = "max_fpb",
  con_thr = FALSE,
  predict_area = NULL
)
#> Predicting individual models

# a list of models
list_p <- sdm_predict(
  models = list(mglm, mraf, mgbm),
  pred = somevar,
  thr = "max_fpb",
  con_thr = FALSE,
  predict_area = NULL
)
#> Predicting list of individual models

# Predict an ensemble model
# (only is possilbe use one fit_ensemble)
ensemble_p <- sdm_predict(
  models = mensemble,
  pred = somevar,
  thr = "max_fpb",
  con_thr = FALSE,
  predict_area = NULL
)
#> Predicting ensembles

# Predict an ensemble of small models
# (only is possible to use one ensemble of small models)
small_p <- sdm_predict(
  models = msmall,
  pred = somevar,
  thr = "max_fpb",
  con_thr = FALSE,
  predict_area = NULL
)
#> Predicting ensemble of small models

## %######################################################%##
#                                                          #
####              Predict model using chunks            ####
#                                                          #
## %######################################################%##
# Predicting models in chunks helps reduce memory requirements in
# cases where models are predicted for large scales and high resolution

ind_p <- sdm_predict(
  models = mglm,
  pred = somevar,
  thr = "max_fpb",
  con_thr = FALSE,
  predict_area = NULL,
  nchunk = 4
)
#> Predicting individual models
# }
```
