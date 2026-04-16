# Fit and validate Domain models

Fit and validate Domain models

## Usage

``` r
fit_dom(
  data,
  response,
  predictors,
  predictors_f = NULL,
  partition = NULL,
  thr = NULL,
  fit_formula = NULL
)
```

## Arguments

- data:

  data.frame. Database with response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1).

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). Usage predictors = c("aet",
  "cwd", "tmin")

- predictors_f:

  character. Vector with the column names of qualitative predictor
  variables (i.e. ordinal or nominal variables; factors). Usage
  predictors_f = c("landform")

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1).
  This is useful for threshold-dependent performance metrics. It is
  possible to use more than one threshold type. It is necessary to
  provide a vector for this argument. The following threshold criteria
  are available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the sensitivity and specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If a sensitivity value is not
    specified, the default used is 0.9.

  If more than one threshold type is used they must be concatenated,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all threshold types
  if none is specified.

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)). Note
  that the variables used here must be consistent with those used in
  response, predictors, and predictors_f arguments

- n_cores:

  numeric. Number of cores to use for parallel processing when metric is
  "domain". Default 1 (no parallelization).

## Value

A list object with:

- model: A tibble with presences. This object can be used for
  predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in the argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- data_ens: Predicted suitability for each test partition. This database
  is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## Details

This function fits and validates Domain models. The Domain model is a
simple model that uses the Gower distance to calculate environmental
similarity between the presence data and test data (Carpenter et al.,
1993). Gower range of values area based on presences data. Gower
distance are transformed to max(0, 1 - Gower). This involves subtracting
the distance from 1 and then ensuring the result is not negative
(clamping it at zero). Gower distance is calculated with
[`map_env_dist`](https://sjevelazco.github.io/flexsdm/reference/map_env_dist.md)
function

This function fit and validate Domain models. The Domain model is a
simple model that uses the Gower distance to calculate the similarity
between the presences training and presence-absences test data.

## References

- Carpenter, G., Gillison, A.N., Winter, J., 1993. DOMAIN: a flexible
  modelling procedure for mapping potential distributions of plants and
  animals. Biodiversity & Conservation 2, 667–680

## See also

[`fit_gam`](https://sjevelazco.github.io/flexsdm/reference/fit_gam.md),
[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md),
[`fit_gbm`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md),
[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md),
[`fit_max`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md),
[`fit_net`](https://sjevelazco.github.io/flexsdm/reference/fit_net.md),
[`fit_raf`](https://sjevelazco.github.io/flexsdm/reference/fit_raf.md),
and
[`fit_svm`](https://sjevelazco.github.io/flexsdm/reference/fit_svm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)

data("spp")
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Extract data
some_sp <- spp %>%
  filter(species == "sp2")

some_sp <-
  sdm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = somevar
  )

# Partition
some_sp <- part_random(
  data = some_sp,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 3, replicates = 5)
)


## %######################################################%##
#                                                          #
####                Fit a Domain model                  ####
#                                                          #
## %######################################################%##
# Fit some models
mdom <- fit_dom(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  predictors_f = NULL,
  fit_formula = NULL,
  partition = ".part",
  thr = c("max_sens_spec"),
  n_cores = 1
)

mdom

# Predict model
ind_p <- sdm_predict(
  models = mdom,
  pred = somevar,
  thr = "max_sens_spec",
  con_thr = TRUE,
  predict_area = NULL
)
plot(ind_p$dom)

## %######################################################%##
#                                                          #
####             Explore Domain suitabiltiy             ####
####             in the environmental space             ####
#                                                          #
## %######################################################%##

p_extra(
  training_data = some_sp %>% dplyr::filter(pr_ab == 1), # select only presences
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = ind_p$dom$dom,
  projection_data = somevar,
  geo_space = FALSE,
  prop_points = 0.3,
  alpha_p = 0.8,
  color_p = "black",
  color_gradient = c("#000033", "#1400FF", "#C729D6", "#FF9C63", "#FFFF60")
)

p_extra(
  training_data = some_sp %>% dplyr::filter(pr_ab == 1), # select only presences
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  predictors = c("CFP_1", "CFP_2"), # Just the first two predictors
  extra_suit_data = ind_p$dom$dom > 0.96, # a binary map
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.4,
  alpha_p = 0.8,
  color_p = "black",
  color_gradient = c("#1400FF", "#C729D6")
)
} # }
```
