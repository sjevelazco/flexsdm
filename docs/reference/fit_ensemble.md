# Ensemble model fitting and validation

Ensemble model fitting and validation

## Usage

``` r
fit_ensemble(
  models,
  ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
  thr = NULL,
  thr_model = NULL,
  metric = NULL
)
```

## Arguments

- models:

  list. A list of models fitted with fit\_ or tune\_ function family.
  Models used for ensemble must have the same presences-absences
  records, partition methods, and threshold types.

- ens_method:

  character. Method used to create ensemble of different models. A
  vector must be provided for this argument. For meansup, meanw or
  pcasup method, it is necessary to provide an evaluation metric and
  threshold in 'metric' and 'thr_model' arguments respectively. By
  default all of the following ensemble methods will be performed:

  - mean: Simple average of the different models.

  - meanw: Weighted average of models based on their performance. An
    evaluation metric and threshold type must be provided.

  - meansup: Average of the best models (those with the evaluation
    metric above the average). An evaluation metric must be provided.

  - meanthr: Averaging performed only with those cells with suitability
    values above the selected threshold.

  - median: Median of the different models.

  Usage ensemble = "meanthr". If several ensemble methods are to be
  implemented it is necessary to concatenate them, e.g., ensemble =
  c("meanw", "meanthr", "median")

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1).
  It is useful for threshold-dependent performance metrics. It is
  possible to use more than one threshold criterion. A vector must be
  provided for this argument. The following threshold criteria are
  available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the sensitivity and specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which Jaccard is the highest.

  - max_sorensen: The threshold at which Sorensen is highest.

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If a sensitivity values is not
    specified, default is 0.9.

  In the case of using more than one threshold type it is necessary
  concatenate threshold types, e.g., thr=c('lpt', 'max_sens_spec',
  'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity',
  sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function
  will use all thresholds if no threshold is specified.

- thr_model:

  character. This threshold is needed for conduct meanw, meandsup, and
  meanthr ensemble methods. It is mandatory to use only one threshold,
  and this must be the same threshold used to fit all the models used in
  the "models" argument. Usage thr_model = 'equal_sens_spec'

- metric:

  character. Performance metric used for selecting the best combination
  of hyper-parameter values. One of the following metrics can be used:
  SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, IMAE, and BOYCE. Default TSS.
  Usage metric = BOYCE

## Value

A list object with:

- models: A list of models used for performing ensemble.

- thr_metric: Threshold and metric specified in the function.

- predictors: A tibble of quantitative (column names with c) and
  qualitative (column names with f) variables used in each models.

- performance: A tibble with performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Those metrics that are threshold-dependent are calculated based on the
  threshold specified in the argument.

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)

# Environmental variables
somevar <-
  system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Species occurrences
data("spp")
set.seed(1)
some_sp <- spp %>%
  dplyr::filter(species == "sp2") %>%
  sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = names(somevar),
    filter_na = TRUE
  ) %>%
  part_random(
    data = .,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )


# gam
mglm <- fit_glm(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  poly = 2
)
mraf <- fit_raf(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
)
mgbm <- fit_gbm(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part"
)

# Fit and validate ensemble model
mensemble <- fit_ensemble(
  models = list(mglm, mraf, mgbm),
  ens_method = "meansup",
  thr = NULL,
  thr_model = "max_sens_spec",
  metric = "TSS"
)

mensemble
} # }
```
