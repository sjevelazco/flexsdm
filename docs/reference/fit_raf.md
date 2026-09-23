# Fit and validate Random Forests models

Fit and validate Random Forests models

## Usage

``` r
fit_raf(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  thr = NULL,
  mtry = sqrt(length(c(predictors, predictors_f))),
  ntree = 500
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
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1),
  needed for threshold-dependent performance metrics. More than one
  threshold type can be used. It is necessary to provide a vector for
  this argument. The following threshold criteria are available:

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
    refers to sensitivity value. If it is not specified a sensitivity
    values, function will use by default 0.9

  If more than one threshold type is used they must be concatenated,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- mtry:

  numeric. Number of variables randomly sampled as candidates at each
  split. Default sqrt(length(c(predictors, predictors_f)))

- ntree:

  numeric. Number of trees to grow. Default 500

## Value

A list object with:

- model: A "randomForest" class object from randomForest package. This
  object can be used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in the argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- data_ens: Predicted suitability for each test partition. This database
  is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## See also

[`fit_gam`](https://sjevelazco.github.io/flexsdm/reference/fit_gam.md),
[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md),
[`fit_gbm`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md),
[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md),
[`fit_max`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md),
[`fit_net`](https://sjevelazco.github.io/flexsdm/reference/fit_net.md),
and
[`fit_svm`](https://sjevelazco.github.io/flexsdm/reference/fit_svm.md).

## Examples

``` r
# \donttest{
data("abies")

# Using k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 10)
)
abies2
#> # A tibble: 1,400 × 14
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
#> # ℹ 3 more variables: depth <dbl>, landform <fct>, .part <int>

rf_t1 <- fit_raf(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  fit_formula = NULL
)
#> Formula used for model fitting:
#> pr_ab ~ aet + ppt_jja + pH + awc + depth + landform
#> Replica number: 1/1
#> Partition number: 1/10
#> Partition number: 2/10
#> Partition number: 3/10
#> Partition number: 4/10
#> Partition number: 5/10
#> Partition number: 6/10
#> Partition number: 7/10
#> Partition number: 8/10
#> Partition number: 9/10
#> Partition number: 10/10

rf_t1$model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = TRUE, ) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 2
#> 
#>         OOB estimate of  error rate: 17.71%
#> Confusion matrix:
#>     0   1 class.error
#> 0 562 138   0.1971429
#> 1 110 590   0.1571429
rf_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
rf_t1$performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 raf   equal_sens_sp…     0.452         700        700    0.817 0.0466    0.817
#> 2 raf   max_sens_spec      0.474         700        700    0.851 0.0926    0.827
#> 3 raf   max_sorensen       0.474         700        700    0.907 0.0443    0.76 
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
rf_t1$performance_part
#> # A tibble: 30 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       1         raf   max_sorensen       0.348          70         70 0.914
#>  2 1       1         raf   max_sens_spec      0.532          70         70 0.786
#>  3 1       1         raf   equal_sens_sp…     0.476          70         70 0.786
#>  4 1       2         raf   max_sorensen       0.492          70         70 0.9  
#>  5 1       2         raf   max_sens_spec      0.492          70         70 0.9  
#>  6 1       2         raf   equal_sens_sp…     0.552          70         70 0.843
#>  7 1       3         raf   max_sorensen       0.472          70         70 0.929
#>  8 1       3         raf   max_sens_spec      0.512          70         70 0.914
#>  9 1       3         raf   equal_sens_sp…     0.574          70         70 0.857
#> 10 1       4         raf   max_sorensen       0.416          70         70 0.957
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
rf_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <fct> <dbl>
#>  1 7      .part      1     0     0.03 
#>  2 26     .part      1     0     0.448
#>  3 40     .part      1     0     0.5  
#>  4 58     .part      1     0     0.06 
#>  5 66     .part      1     0     0.182
#>  6 67     .part      1     0     0.474
#>  7 74     .part      1     0     0.152
#>  8 81     .part      1     0     0.078
#>  9 87     .part      1     0     0.868
#> 10 94     .part      1     0     0.644
#> # ℹ 1,390 more rows

# Using bootstrap partition method and only with presence-absence
# and get performance for several method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 10, proportion = 0.7)
)
abies2
#> # A tibble: 1,400 × 23
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
#> # ℹ 12 more variables: depth <dbl>, landform <fct>, .part1 <chr>, .part2 <chr>,
#> #   .part3 <chr>, .part4 <chr>, .part5 <chr>, .part6 <chr>, .part7 <chr>,
#> #   .part8 <chr>, .part9 <chr>, .part10 <chr>

rf_t2 <- fit_raf(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  fit_formula = NULL,
  mtry = 2,
  ntree = 500
)
#> Formula used for model fitting:
#> pr_ab ~ aet + ppt_jja + pH + awc + depth + landform
#> Replica number: 1/10
#> Partition number: 1/1
#> Replica number: 2/10
#> Partition number: 1/1
#> Replica number: 3/10
#> Partition number: 1/1
#> Replica number: 4/10
#> Partition number: 1/1
#> Replica number: 5/10
#> Partition number: 1/1
#> Replica number: 6/10
#> Partition number: 1/1
#> Replica number: 7/10
#> Partition number: 1/1
#> Replica number: 8/10
#> Partition number: 1/1
#> Replica number: 9/10
#> Partition number: 1/1
#> Replica number: 10/10
#> Partition number: 1/1
rf_t2
#> $model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = TRUE, ) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 2
#> 
#>         OOB estimate of  error rate: 17.71%
#> Confusion matrix:
#>     0   1 class.error
#> 0 562 138   0.1971429
#> 1 110 590   0.1571429
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
#> 
#> $performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 raf   equal_sens_sp…     0.452         700        700    0.931 0.0117    0.931
#> 2 raf   max_sens_spec      0.474         700        700    0.950 0.0183    0.929
#> 3 raf   max_sorensen       0.474         700        700    0.958 0.0212    0.92 
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 30 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       1         raf   max_sorensen       0.374         210        210 0.981
#>  2 1       1         raf   max_sens_spec      0.444         210        210 0.948
#>  3 1       1         raf   equal_sens_sp…     0.51          210        210 0.924
#>  4 2       1         raf   max_sorensen       0.58          210        210 0.933
#>  5 2       1         raf   max_sens_spec      0.59          210        210 0.929
#>  6 2       1         raf   equal_sens_sp…     0.538         210        210 0.938
#>  7 3       1         raf   max_sorensen       0.384         210        210 0.967
#>  8 3       1         raf   max_sens_spec      0.384         210        210 0.967
#>  9 3       1         raf   equal_sens_sp…     0.472         210        210 0.910
#> 10 4       1         raf   max_sorensen       0.448         210        210 0.967
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
#> $data_ens
#> # A tibble: 4,200 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <fct> <dbl>
#>  1 1      .part1     1     0     0.286
#>  2 2      .part1     1     0     0.028
#>  3 6      .part1     1     0     0.208
#>  4 13     .part1     1     0     0.108
#>  5 16     .part1     1     0     0.028
#>  6 17     .part1     1     0     0.138
#>  7 18     .part1     1     0     0    
#>  8 20     .part1     1     0     0.38 
#>  9 24     .part1     1     0     0.316
#> 10 25     .part1     1     0     0.122
#> # ℹ 4,190 more rows
#> 
# }
```
