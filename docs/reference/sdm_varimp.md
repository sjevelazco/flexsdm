# Calculate permutation-based variable importance scores for SDMs

This function calculates variable importance scores for species
distribution models (SDMs) based on permutation-based approach.

## Usage

``` r
sdm_varimp(
  models,
  data,
  response,
  predictors,
  n_sim = 50,
  n_cores = 1,
  thr = NULL,
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

- data:

  data.frame. Database with response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1).

- predictors:

  character. Vector with the column names of predictor variables. Usage
  predictors = c("aet", "cwd", "tmin")

- n_sim:

  integer. The number of Monte Carlo replications to perform. Default
  is 50. The results from each replication are averaged together (the
  standard deviation will also be returned).

- n_cores:

  numeric. Number of cores use for parallelization. Default 1

- thr:

  character. Threshold criterion used to get binary suitability values
  (i.e. 0,1). Used for threshold-dependent performance metrics. It is
  possible to use more than one threshold type. A vector must be
  provided for this argument. The following threshold criteria are
  available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the Sensitivity and Specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the Sensitivity and
    Specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is the
    highest.

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is the highest.

  - sensitivity: Threshold based on a specified Sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to Sensitivity value. If a sensitivity value is not
    specified, the default value is 0.9

  \#' If more than one threshold type is used, concatenate threshold
  types, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or
  thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or
  thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all
  thresholds if no threshold type is specified

- clamp:

  logical. If TRUE, predictors and features are restricted to the range
  seen during model training.

- pred_type:

  character. Type of response required available "link", "exponential",
  "cloglog" and "logistic". Default "cloglog"

## Value

a tibble with the columns:

- model: model name

- threshold: threshold names

- predictors: predictor names

- from TPR to IMAE: performance metrics

## Details

This function calculates variable importance scores for species
distribution models (SDMs) based on permutation-based approach. Thus,
the function calculates the model performance using the original data
and the permuted data. The difference between the two performances is
the variable importance score.

## Examples

``` r
if (FALSE) { # \dontrun{
require(tidyr)
require(dplyr)

data(abies)
abies

data(backg)
backg

# In this example we will partition the data using the k-fold method

abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

max_t1 <- fit_max(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  background = backg2,
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  clamp = TRUE,
  classes = "default",
  pred_type = "cloglog",
  regmult = 1
)

net_t1 <- fit_net(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
)

svm_f1 <- fit_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
)

vip_t <- sdm_varimp(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "landform"),
  models = list(max_t1, net_t1, svm_f1),
  clamp = TRUE,
  pred_type = "cloglog",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  n_sim = 50,
  n_cores = 5
)

vip_t

# Plot the variable importance for AUC TSS and SORENSEN for the
# threshold that maximizes Sorensen metric and Maxent
vip_t %>%
  pivot_longer(
    cols = TPR:IMAE,
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::filter(threshold == "max_sorensen") %>%
  dplyr::filter(metric %in% c("AUC", "TSS", "SORENSEN")) %>%
  dplyr::filter(model == "max") %>%
  ggplot(aes(x = reorder(predictors, value), y = value, fill = predictors)) +
  geom_col(
    col = "black",
    show.legend = FALSE
  ) +
  facet_wrap(~metric, scales = "free_x") +
  labs(x = "Predictors", y = "Variable Importance") +
  theme_classic() +
  coord_flip()
} # }
```
