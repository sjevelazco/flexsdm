# Fit and validate Neural Networks models with exploration of hyper-parameters

Fit and validate Neural Networks models with exploration of
hyper-parameters

## Usage

``` r
tune_net(
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
  that the variable names used here must be consistent with those used
  in response, predictors, and predictors_f arguments. Defaul NULL.

- partition:

  character. Column name with training and validation partition groups.

- grid:

  data.frame. A data frame object with algorithm hyper-parameters values
  to be tested. It is recommended to generate this data.frame with the
  grid() function.

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

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If it is not specified a sensitivity
    values, function will use by default 0.9.

  If using more than one threshold type concatenate them, e.g.,
  thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- metric:

  character. Performance metric used for selecting the best combination
  of hyper-parameter values. One of the following metrics can be used:
  SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as
  default.

- n_cores:

  numeric. Number of cores use for parallelization. Default 1

## Value

A list object with:

- model: A "nnet" class object from nnet package. This object can be
  used for predicting.

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
[`tune_raf`](https://sjevelazco.github.io/flexsdm/reference/tune_raf.md),
and
[`tune_svm`](https://sjevelazco.github.io/flexsdm/reference/tune_svm.md).

## Examples

``` r
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

# Partitioning the data with the k-fold method

abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

# pr_ab columns is species presence and absences (i.e. the response variable)
# from aet to landform are the predictors variables (landform is a qualitative variable)

# Hyper-parameter values for tuning
tune_grid <-
  expand.grid(
    size = c(2, 4, 6, 8, 10),
    decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10)
  )

net_t <-
  tune_net(
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
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5

# Outputs
net_t$model
#> a 22-8-1 network with 193 weights
#> inputs: aet cwd tmin ppt_djf ppt_jja pH awc depth landform2 landform3 landform4 landform5 landform6 landform7 landform8 landform9 landform10 landform11 landform12 landform13 landform14 landform15 
#> output(s): pr_ab 
#> options were - entropy fitting  decay=1
net_t$predictors
#> # A tibble: 1 × 9
#>   c1    c2    c3    c4      c5      c6    c7    c8    f       
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth landform
net_t$performance
#> # A tibble: 1 × 35
#>    size decay model threshold   thr_value n_presences n_absences TPR_mean TPR_sd
#>   <dbl> <dbl> <chr> <chr>           <dbl>       <int>      <int>    <dbl>  <dbl>
#> 1     8     1 net   max_sens_s…     0.457         700        700    0.904 0.0297
#> # ℹ 26 more variables: TNR_mean <dbl>, TNR_sd <dbl>, W_TPR_TNR_mean <dbl>,
#> #   W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,
#> #   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,
#> #   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>,
#> #   KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
net_t$performance_part
#> # A tibble: 5 × 21
#>   replica partition model threshold thr_value n_presences n_absences   TPR   TNR
#>   <chr>   <chr>     <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#> 1 1       1         net   max_sens…     0.590         140        140 0.893 0.886
#> 2 1       2         net   max_sens…     0.385         140        140 0.943 0.864
#> 3 1       3         net   max_sens…     0.535         140        140 0.929 0.871
#> 4 1       4         net   max_sens…     0.643         140        140 0.879 0.864
#> 5 1       5         net   max_sens…     0.645         140        140 0.879 0.807
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
net_t$hyper_performance
#> # A tibble: 40 × 32
#>     size  decay model threshold   TPR_mean TPR_sd TNR_mean TNR_sd W_TPR_TNR_mean
#>    <dbl>  <dbl> <chr> <chr>          <dbl>  <dbl>    <dbl>  <dbl>          <dbl>
#>  1     2  0.001 net   max_sens_s…    1     0         0     0               0.5  
#>  2     2  0.05  net   max_sens_s…    0.881 0.0854    0.779 0.124           0.83 
#>  3     2  0.1   net   max_sens_s…    0.883 0.0697    0.833 0.0499          0.858
#>  4     2  1     net   max_sens_s…    0.887 0.0561    0.81  0.0465          0.849
#>  5     2  3     net   max_sens_s…    0.893 0.0368    0.826 0.0439          0.859
#>  6     2  4     net   max_sens_s…    0.84  0.0894    0.811 0.0963          0.826
#>  7     2  5     net   max_sens_s…    0.913 0.0653    0.731 0.112           0.822
#>  8     2 10     net   max_sens_s…    0.919 0.0383    0.736 0.0991          0.827
#>  9     4  0.001 net   max_sens_s…    0.804 0.132     0.693 0.236           0.749
#> 10     4  0.05  net   max_sens_s…    0.883 0.0318    0.834 0.0348          0.859
#> # ℹ 30 more rows
#> # ℹ 23 more variables: W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>,
#> #   SORENSEN_sd <dbl>, JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>,
#> #   FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>,
#> #   KAPPA_mean <dbl>, KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>,
#> #   AUC_mean <dbl>, AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>,
#> #   CRPS_mean <dbl>, CRPS_sd <dbl>, IMAE_mean <dbl>, IMAE_sd <dbl>
net_t$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab    pred
#>    <chr>  <chr>      <chr> <fct>   <dbl>
#>  1 3      .part      1     0     0.168  
#>  2 5      .part      1     0     0.393  
#>  3 7      .part      1     0     0.0542 
#>  4 20     .part      1     0     0.209  
#>  5 22     .part      1     0     0.0881 
#>  6 26     .part      1     0     0.756  
#>  7 29     .part      1     0     0.00417
#>  8 30     .part      1     0     0.00420
#>  9 51     .part      1     0     0.0783 
#> 10 62     .part      1     0     0.209  
#> # ℹ 1,390 more rows
# }
```
