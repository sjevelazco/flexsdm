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
# \donttest{
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
#> 6 rows were excluded from database because NAs were found


# gam
mglm <- fit_glm(
  data = some_sp,
  response = "pr_ab",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  poly = 2
)
#> Formula used for model fitting:
#> pr_ab ~ CFP_1 + CFP_2 + CFP_3 + CFP_4 + I(CFP_1^2) + I(CFP_2^2) + I(CFP_3^2) + I(CFP_4^2)
#> Replica number: 1/1
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
#> Replica number: 1/1
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
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

# Fit and validate ensemble model
mensemble <- fit_ensemble(
  models = list(mglm, mraf, mgbm),
  ens_method = "meansup",
  thr = NULL,
  thr_model = "max_sens_spec",
  metric = "TSS"
)
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |======================================================================| 100%

mensemble
#> $models
#> $models$m_1
#> $models$m_1$model
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)        CFP_1        CFP_2        CFP_3        CFP_4   I(CFP_1^2)  
#>   2.307e+01   -6.272e-02    9.755e-01    4.675e-02    9.784e-01    3.183e-05  
#>  I(CFP_2^2)   I(CFP_3^2)   I(CFP_4^2)  
#>  -3.904e-02   -4.391e-04   -9.860e-02  
#> 
#> Degrees of Freedom: 93 Total (i.e. Null);  85 Residual
#> Null Deviance:       108.9 
#> Residual Deviance: 56.05     AIC: 74.05
#> 
#> $models$m_1$performance
#> # A tibble: 7 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 glm   equal_sens_sp…    0.413           25         69    0.801 0.0656    0.812
#> 2 glm   lpt               0.0707          25         69    1     0         0.681
#> 3 glm   max_fpb           0.554           25         69    0.764 0.206     0.942
#> 4 glm   max_jaccard       0.554           25         69    0.764 0.206     0.942
#> 5 glm   max_sens_spec     0.459           25         69    0.963 0.0642    0.783
#> 6 glm   max_sorensen      0.554           25         69    0.764 0.206     0.942
#> 7 glm   sensitivity       0.246           25         69    1     0         0.681
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> 
#> $models$m_2
#> $models$m_2$model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = TRUE, ) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 2
#> 
#>         OOB estimate of  error rate: 11.7%
#> Confusion matrix:
#>    0  1 class.error
#> 0 63  6  0.08695652
#> 1  5 20  0.20000000
#> 
#> $models$m_2$performance
#> # A tibble: 7 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 raf   equal_sens_sp…     0.684          25         69    0.843 0.0561    0.826
#> 2 raf   lpt                0.684          25         69    1     0         0.725
#> 3 raf   max_fpb            0.684          25         69    0.806 0.120     0.928
#> 4 raf   max_jaccard        0.684          25         69    0.806 0.120     0.928
#> 5 raf   max_sens_spec      0.684          25         69    0.884 0.111     0.870
#> 6 raf   max_sorensen       0.684          25         69    0.806 0.120     0.928
#> 7 raf   sensitivity        0.698          25         69    1     0         0.725
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> 
#> $models$m_3
#> $models$m_3$model
#> gbm::gbm(formula = formula1, distribution = "bernoulli", data = data, 
#>     n.trees = n_trees, n.minobsinnode = n_minobsinnode, shrinkage = shrinkage)
#> A gradient boosted model with bernoulli loss function.
#> 100 iterations were performed.
#> There were 4 predictors of which 4 had non-zero influence.
#> 
#> $models$m_3$performance
#> # A tibble: 7 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 gbm   equal_sens_sp…     0.600          25         69    0.843 0.0561    0.855
#> 2 gbm   lpt                0.232          25         69    1     0         0.826
#> 3 gbm   max_fpb            0.671          25         69    0.926 0.128     0.899
#> 4 gbm   max_jaccard        0.671          25         69    0.926 0.128     0.899
#> 5 gbm   max_sens_spec      0.232          25         69    1     0         0.826
#> 6 gbm   max_sorensen       0.671          25         69    0.926 0.128     0.899
#> 7 gbm   sensitivity        0.319          25         69    1     0         0.826
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> 
#> 
#> $thr_metric
#> [1] "max_sens_spec" "TSS_mean"     
#> 
#> $predictors
#> # A tibble: 3 × 4
#>   c1    c2    c3    c4   
#>   <chr> <chr> <chr> <chr>
#> 1 CFP_1 CFP_2 CFP_3 CFP_4
#> 2 CFP_1 CFP_2 CFP_3 CFP_4
#> 3 CFP_1 CFP_2 CFP_3 CFP_4
#> 
#> $performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 meansup equal_sens_…     0.347          25         69    0.843 0.0561    0.855
#> 2 meansup lpt              0.131          25         69    1     0         0.826
#> 3 meansup max_fpb          0.518          25         69    0.926 0.128     0.899
#> 4 meansup max_jaccard      0.518          25         69    0.926 0.128     0.899
#> 5 meansup max_sens_sp…     0.518          25         69    1     0         0.826
#> 6 meansup max_sorensen     0.518          25         69    0.926 0.128     0.899
#> 7 meansup sensitivity      0.271          25         69    1     0         0.826
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 21 × 21
#>    model   replica partition threshold    thr_value n_presences n_absences   TPR
#>    <chr>   <chr>   <chr>     <chr>            <dbl>       <int>      <int> <dbl>
#>  1 meansup 1       1         max_sorensen     0.518           9         23 0.778
#>  2 meansup 1       1         max_jaccard      0.518           9         23 0.778
#>  3 meansup 1       1         max_fpb          0.518           9         23 0.778
#>  4 meansup 1       1         max_sens_sp…     0.186           9         23 1    
#>  5 meansup 1       1         equal_sens_…     0.314           9         23 0.778
#>  6 meansup 1       1         lpt              0.186           9         23 1    
#>  7 meansup 1       1         sensitivity      0.186           9         23 1    
#>  8 meansup 1       2         max_sorensen     0.131           8         23 1    
#>  9 meansup 1       2         max_jaccard      0.131           8         23 1    
#> 10 meansup 1       2         max_fpb          0.131           8         23 1    
#> # ℹ 11 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
# }
```
