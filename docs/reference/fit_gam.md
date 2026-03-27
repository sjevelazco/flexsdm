# Fit and validate Generalized Additive Models

Fit and validate Generalized Additive Models

## Usage

``` r
fit_gam(
  data,
  response,
  predictors,
  predictors_f = NULL,
  select_pred = FALSE,
  partition = NULL,
  thr = NULL,
  fit_formula = NULL,
  k = -1
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

- select_pred:

  logical. Perform predictor selection. Default FALSE.

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

- k:

  integer. The dimension of the basis used to represent the smooth term.
  Default -1 (i.e., k=10). See the help in ?mgcv::s.

## Value

A list object with:

- model: A "gam" class object from mgcv package. This object can be used
  for predicting.

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

This function fits GAM using mgvc package, with Binomial distribution
family and thin plate regression spline as a smoothing basis (see
?mgvc::s).

## See also

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
data("abies")

# Using k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 10)
)
abies2

gam_t1 <- fit_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  select_pred = FALSE,
  partition = ".part",
  thr = "max_sens_spec"
)
gam_t1$model
gam_t1$predictors
gam_t1$performance
gam_t1$performance_part

# Specifying the formula explicitly
require(mgcv)
gam_t2 <- fit_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  select_pred = FALSE,
  partition = ".part",
  thr = "max_sens_spec",
  fit_formula = stats::formula(pr_ab ~ s(aet) +
    s(ppt_jja) +
    s(pH) + landform)
)

gam_t2$model
gam_t2$predictors
gam_t2$performance %>% dplyr::select(ends_with("_mean"))

# Using repeated k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 5)
)
abies2

gam_t3 <- fit_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("ppt_jja", "pH", "awc"),
  predictors_f = c("landform"),
  select_pred = FALSE,
  partition = ".part",
  thr = "max_sens_spec"
)
gam_t3
} # }
```
