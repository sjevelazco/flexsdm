# Fit and validate Random Forest models with exploration of hyper-parameters that optimize performance

Fit and validate Random Forest models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_raf(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  grid = NULL,
  thr = NULL,
  metric = "TSS",
  n_cores = 1
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
  variables (i.e. ordinal or nominal variables type). Usage predictors_f
  = c("landform")

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)). Note
  that the variables used here must be consistent with those used in
  response, predictors, and predictors_f arguments. Default NULL

- partition:

  character. Column name with training and validation partition groups.

- grid:

  data.frame. A data frame object with algorithm hyper-parameters values
  to be tested. It is recommended to generate this data.frame with the
  grid() function. Hyper-parameter needed for tuning is 'mtry' and
  'ntree'. The maximum mtry cannot exceed the total number of
  predictors.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1),
  needed for threshold-dependent performance metrics. It is possible to
  use more than one threshold type. It is necessary to provide a vector
  for this argument. The following threshold types are available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the sensitivity and specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which FPB is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If it is not specified a sensitivity
    values, function will use by default 0.9

  If using more than one threshold type concatenate them, e.g.,
  thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- metric:

  character. Performance metric used for selecting the best combination
  of hyper -parameter values. One of the following metrics can be used:
  SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as
  default.

- n_cores:

  numeric. Number of cores use for parallelization. Default 1

## Value

A list object with:

- model: A "randomForest" class object from randomForest package. This
  object can be used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Hyper-parameters values and performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md))
  for the best hyper-parameters combination.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- hyper_performance: Performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md))
  for each combination of the hyper-parameters.

- data_ens: Predicted suitability for each test partition based on the
  best model. This database is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## See also

[`tune_gbm`](https://sjevelazco.github.io/flexsdm/reference/tune_gbm.md),
[`tune_max`](https://sjevelazco.github.io/flexsdm/reference/tune_max.md),
[`tune_net`](https://sjevelazco.github.io/flexsdm/reference/tune_net.md),
and
[`tune_svm`](https://sjevelazco.github.io/flexsdm/reference/tune_svm.md).

## Examples

``` r
# Toy data: a small subset of 'abies' to keep this example fast
data("abies")
abies_toy <- abies[c(1:30, 701:730), ]
abies_toy <- part_random(
  data = abies_toy,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

tune_toy <- tune_raf(
  data = abies_toy,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin"),
  partition = ".part",
  grid = expand.grid(mtry = c(1, 2), ntree = c(50, 100)),
  thr = "max_sens_spec",
  metric = "TSS",
  n_cores = 1
)
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmin
#> Tuning model...
#> Replica number: 1/1
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmin
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
tune_toy$performance
#> # A tibble: 1 × 35
#>    mtry ntree model threshold   thr_value n_presences n_absences TPR_mean TPR_sd
#>   <dbl> <dbl> <chr> <chr>           <dbl>       <int>      <int>    <dbl>  <dbl>
#> 1     1    50 raf   max_sens_s…      0.72          30         30    0.967 0.0577
#> # ℹ 26 more variables: TNR_mean <dbl>, TNR_sd <dbl>, W_TPR_TNR_mean <dbl>,
#> #   W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,
#> #   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,
#> #   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>,
#> #   KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>

# \donttest{
data(abies)
abies
#> # A tibble: 1,400 × 13
#>       id pr_ab        x        y   aet   cwd  tmin ppt_djf ppt_jja    pH    awc
#>    <int> <dbl>    <dbl>    <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl>  <dbl>
#>  1   715     0  -95417.  314240.  323.  546.  1.24    62.7   17.8   5.77 0.108 
#>  2  5680     0   98987. -159415.  448.  815.  9.43   130.     6.43  5.60 0.160 
#>  3  7907     0  121474.  -99463.  182.  271. -4.95   151.    11.2   0    0     
#>  4  1850     0  -39976.  -17456.  372.  946.  8.78   116.     2.70  6.41 0.0972
#>  5  1702     0  111372.  -91404.  209.  399. -4.03   165.     9.27  0    0     
#>  6 10036     0 -255715.  392229.  308.  535.  4.66   166.    16.5   5.70 0.0777
#>  7 12384     0 -311765.  380213.  568.  352.  4.38   480.    41.2   5.80 0.110 
#>  8  6513     0  111360. -120229.  327.  633.  4.93   163.     8.91  1.18 0.0116
#>  9  9884     0 -284326.  442136.  377.  446.  3.99   296.    16.8   5.96 0.0900
#> 10  8651     0  137640. -110538.  215.  265. -4.62   180.     9.57  0    0     
#> # ℹ 1,390 more rows
#> # ℹ 2 more variables: depth <dbl>, landform <fct>

# Partition the data with the k-fold method

abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

tune_grid <- expand.grid(
  mtry = c(2, 4),
  ntree = c(100, 300)
)

tune_grid
#>   mtry ntree
#> 1    2   100
#> 2    4   100
#> 3    2   300
#> 4    4   300

rf_t <-
  tune_raf(
    data = abies2,
    response = "pr_ab",
    predictors = c(
      "aet", "cwd", "tmin", "ppt_djf",
      "ppt_jja", "pH", "awc", "depth"
    ),
    predictors_f = c("landform"),
    partition = ".part",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 1
  )
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmin + ppt_djf + ppt_jja + pH + awc + depth + landform
#> Tuning model...
#> Replica number: 1/1
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmin + ppt_djf + ppt_jja + pH + awc + depth + landform
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

# Outputs
rf_t$model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = TRUE, ) 
#>                Type of random forest: classification
#>                      Number of trees: 300
#> No. of variables tried at each split: 2
#> 
#>         OOB estimate of  error rate: 11.21%
#> Confusion matrix:
#>     0   1 class.error
#> 0 602  98  0.14000000
#> 1  59 641  0.08428571
rf_t$predictors
#> # A tibble: 1 × 9
#>   c1    c2    c3    c4      c5      c6    c7    c8    f       
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth landform
rf_t$performance_part
#> # A tibble: 3 × 21
#>   replica partition model threshold thr_value n_presences n_absences   TPR   TNR
#>   <chr>   <chr>     <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#> 1 1       1         raf   max_sens…     0.53          234        234 0.902 0.876
#> 2 1       2         raf   max_sens…     0.647         233        233 0.906 0.906
#> 3 1       3         raf   max_sens…     0.613         233        233 0.858 0.888
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
rf_t$hyper_performance
#> # A tibble: 4 × 32
#>    mtry ntree model threshold    TPR_mean  TPR_sd TNR_mean TNR_sd W_TPR_TNR_mean
#>   <dbl> <dbl> <chr> <chr>           <dbl>   <dbl>    <dbl>  <dbl>          <dbl>
#> 1     2   100 raf   max_sens_sp…    0.893 0.0188     0.880 0.0113          0.886
#> 2     2   300 raf   max_sens_sp…    0.887 0.0179     0.890 0.0220          0.889
#> 3     4   100 raf   max_sens_sp…    0.900 0.0126     0.871 0.0334          0.886
#> 4     4   300 raf   max_sens_sp…    0.897 0.00430    0.877 0.0322          0.887
#> # ℹ 23 more variables: W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>,
#> #   SORENSEN_sd <dbl>, JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>,
#> #   FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>,
#> #   KAPPA_mean <dbl>, KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>,
#> #   AUC_mean <dbl>, AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>,
#> #   CRPS_mean <dbl>, CRPS_sd <dbl>, IMAE_mean <dbl>, IMAE_sd <dbl>
rf_t$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab    pred
#>    <chr>  <chr>      <chr> <fct>   <dbl>
#>  1 2      .part      1     0     0.0167 
#>  2 4      .part      1     0     0.00333
#>  3 8      .part      1     0     0.107  
#>  4 10     .part      1     0     0.18   
#>  5 15     .part      1     0     0.0233 
#>  6 16     .part      1     0     0.0467 
#>  7 24     .part      1     0     0.9    
#>  8 27     .part      1     0     0.22   
#>  9 29     .part      1     0     0.00333
#> 10 30     .part      1     0     0.0467 
#> # ℹ 1,390 more rows
# }
```
