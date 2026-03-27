# Fit and validate Generalized Boosted Regression models based on Ensembles of Small Models approach

This function constructs Generalized Boosted Regression using the
Ensembles of Small Models (ESM) approach (Breiner et al., 2015, 2018).

## Usage

``` r
esm_gbm(
  data,
  response,
  predictors,
  partition,
  thr = NULL,
  n_trees = 100,
  n_minobsinnode = NULL,
  shrinkage = 0.1
)
```

## Arguments

- data:

  data.frame. Database with the response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1)

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). This can only construct models
  with continuous variables and does not allow categorical variables.
  Usage predictors = c("aet", "cwd", "tmin")

- partition:

  character. Column name with training and validation partition groups.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1).
  It is useful for threshold-dependent performance metrics. It is
  possible to use more than one threshold type. It is necessary to
  provide a vector for this argument. The following threshold criteria
  are available:

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
    specified, the default value is 0.9.

  In the case of use more than one threshold type it is necessary
  concatenate threshold types, e.g., thr=c('max_sens_spec',
  'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'),
  or thr=c('max_sens_spec', 'sensitivity'). Function will use all
  thresholds if no threshold is specified.

- n_trees:

  Integer specifying the total number of trees to fit. This is
  equivalent to the number of iterations and the number of basis
  functions in the additive expansion. Default is 100.

- n_minobsinnode:

  Integer specifying the minimum number of observations in the terminal
  nodes of the trees. Note that this is the actual number of
  observations, not the total weight. If n_minobsinnode is NULL, this
  parameter will assume a value equal to nrow(data)\*0.5/4. Default is
  NULL.

- shrinkage:

  Numeric. This parameter applied to each tree in the expansion. Also
  known as the learning rate or step-size reduction; 0.001 to 0.1
  usually works, but a smaller learning rate typically requires more
  trees. Default is 0.1.

## Value

A list object with:

- esm_model: A list with "gbm" class object from gbm package for each
  bivariate model. This object can be used for predicting ensembles of
  small models with
  [`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)
  function.

- predictors: A tibble with variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in thr argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

## Details

This method consists of creating bivariate models with all the pair-wise
combinations of predictors and perform an ensemble based on the average
of suitability weighted by Somers' D metric (D = 2 x (AUC -0.5)). ESM is
recommended for modeling species with few occurrences. This function
does not allow categorical variables because the use of these types of
variables could be problematic when using with few occurrences. For
further detail see Breiner et al. (2015, 2018).

## References

- Breiner, F. T., Guisan, A., Bergamini, A., & Nobis, M. P. (2015).
  Overcoming limitations of modelling rare species by using ensembles of
  small models. Methods in Ecology and Evolution, 6(10), 1210-218.
  https://doi.org/10.1111/2041-210X.12403

- Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A. (2018).
  Optimizing ensembles of small models for predicting the distribution
  of species with few occurrences. Methods in Ecology and Evolution,
  9(4), 802-808. https://doi.org/10.1111/2041-210X.12957

## See also

[`esm_gam`](https://sjevelazco.github.io/flexsdm/reference/esm_gam.md),
[`esm_gau`](https://sjevelazco.github.io/flexsdm/reference/esm_gau.md),
[`esm_glm`](https://sjevelazco.github.io/flexsdm/reference/esm_glm.md),
[`esm_max`](https://sjevelazco.github.io/flexsdm/reference/esm_max.md),
[`esm_net`](https://sjevelazco.github.io/flexsdm/reference/esm_net.md),
and
[`esm_svm`](https://sjevelazco.github.io/flexsdm/reference/esm_svm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data("abies")
require(dplyr)

# Using k-fold partition method
set.seed(10)
abies2 <- abies %>%
  na.omit() %>%
  group_by(pr_ab) %>%
  dplyr::slice_sample(n = 10) %>%
  group_by()

abies2 <- part_random(
  data = abies2,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 3, replicates = 5)
)
abies2

# Without threshold specification and with kfold
esm_gbm_t1 <- esm_gbm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL,
  n_trees = 100,
  n_minobsinnode = NULL,
  shrinkage = 0.1
)

esm_gbm_t1$esm_model # bivariate model
esm_gbm_t1$predictors
esm_gbm_t1$performance
esm_gbm_t1$performance_part
} # }
```
