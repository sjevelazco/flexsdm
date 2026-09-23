# Fit and validate Generalized Boosted Regression models

Fit and validate Generalized Boosted Regression models

## Usage

``` r
fit_gbm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  thr = NULL,
  n_trees = 100,
  n_minobsinnode = as.integer(nrow(data) * 0.9/4),
  shrinkage = 0.1
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
  response, predictors, and predictors_f arguments. Default is NULL.

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1)
  needed for threshold-dependent performance metrics. It is possible to
  use more than one threshold type. It is necessary to provide a vector
  for this argument. The following threshold criteria are available:

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
    specified, the default used is 0.9

  If more than one threshold type is used they must be concatenated,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- n_trees:

  Integer specifying the total number of trees to fit. This is
  equivalent to the number of iterations and the number of basis
  functions in the additive expansion. Default is 100.

- n_minobsinnode:

  Integer specifying the minimum number of observations in the terminal
  nodes of the trees. Note that this is the actual number of
  observations, not the total weight. The default value used is
  nrow(data)\*0.5/4

- shrinkage:

  Numeric. This parameter applied to each tree in the expansion. Also
  known as the learning rate or step-size reduction; 0.001 to 0.1
  usually works, but a smaller learning rate typically requires more
  trees. Default is 0.1.

## Value

A list object with:

- model: A "gbm" class object from gbm package. This object can be used
  for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in thr argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- data_ens: Predicted suitability for each test partition based on the
  best model. This database is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## See also

[`fit_gam`](https://sjevelazco.github.io/flexsdm/reference/fit_gam.md),
[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md),
[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md),
[`fit_max`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md),
[`fit_net`](https://sjevelazco.github.io/flexsdm/reference/fit_net.md),
[`fit_raf`](https://sjevelazco.github.io/flexsdm/reference/fit_raf.md),
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

gbm_t1 <- fit_gbm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
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
gbm_t1$model
#> gbm::gbm(formula = formula1, distribution = "bernoulli", data = data, 
#>     n.trees = n_trees, n.minobsinnode = n_minobsinnode, shrinkage = shrinkage)
#> A gradient boosted model with bernoulli loss function.
#> 100 iterations were performed.
#> There were 6 predictors of which 6 had non-zero influence.
gbm_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
gbm_t1$performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 gbm   equal_sens_sp…     0.508         700        700    0.681 0.0415    0.684
#> 2 gbm   max_sens_spec      0.473         700        700    0.757 0.125     0.671
#> 3 gbm   max_sorensen       0.394         700        700    0.904 0.0500    0.484
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
gbm_t1$performance_part
#> # A tibble: 30 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       1         gbm   max_sorensen       0.430          70         70 0.814
#>  2 1       1         gbm   max_sens_spec      0.436          70         70 0.786
#>  3 1       1         gbm   equal_sens_sp…     0.513          70         70 0.643
#>  4 1       2         gbm   max_sorensen       0.457          70         70 0.871
#>  5 1       2         gbm   max_sens_spec      0.475          70         70 0.829
#>  6 1       2         gbm   equal_sens_sp…     0.586          70         70 0.6  
#>  7 1       3         gbm   max_sorensen       0.391          70         70 0.943
#>  8 1       3         gbm   max_sens_spec      0.391          70         70 0.943
#>  9 1       3         gbm   equal_sens_sp…     0.491          70         70 0.7  
#> 10 1       4         gbm   max_sorensen       0.369          70         70 0.957
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
gbm_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <dbl> <dbl>
#>  1 4      .part      1         0 0.151
#>  2 13     .part      1         0 0.595
#>  3 19     .part      1         0 0.388
#>  4 24     .part      1         0 0.776
#>  5 31     .part      1         0 0.403
#>  6 48     .part      1         0 0.358
#>  7 51     .part      1         0 0.595
#>  8 62     .part      1         0 0.213
#>  9 79     .part      1         0 0.336
#> 10 87     .part      1         0 0.845
#> # ℹ 1,390 more rows

# Using bootstrap partition method
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

gbm_t2 <- fit_gbm(
  data = abies2,
  response = "pr_ab",
  predictors = c("ppt_jja", "pH", "awc"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = "max_sens_spec"
)
#> Formula used for model fitting:
#> pr_ab ~ ppt_jja + pH + awc + landform
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
gbm_t2
#> $model
#> gbm::gbm(formula = formula1, distribution = "bernoulli", data = data, 
#>     n.trees = n_trees, n.minobsinnode = n_minobsinnode, shrinkage = shrinkage)
#> A gradient boosted model with bernoulli loss function.
#> 100 iterations were performed.
#> There were 4 predictors of which 4 had non-zero influence.
#> 
#> $predictors
#> # A tibble: 1 × 4
#>   c1      c2    c3    f       
#>   <chr>   <chr> <chr> <chr>   
#> 1 ppt_jja pH    awc   landform
#> 
#> $performance
#> # A tibble: 1 × 33
#>   model threshold     thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 gbm   max_sens_spec     0.535         700        700    0.668 0.0590    0.706
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 10 × 21
#>    replica partition model threshold     thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>             <dbl>       <int>      <int> <dbl>
#>  1 1       1         gbm   max_sens_spec     0.524         210        210 0.629
#>  2 2       1         gbm   max_sens_spec     0.471         210        210 0.7  
#>  3 3       1         gbm   max_sens_spec     0.560         210        210 0.605
#>  4 4       1         gbm   max_sens_spec     0.490         210        210 0.710
#>  5 5       1         gbm   max_sens_spec     0.471         210        210 0.748
#>  6 6       1         gbm   max_sens_spec     0.500         210        210 0.652
#>  7 7       1         gbm   max_sens_spec     0.536         210        210 0.614
#>  8 8       1         gbm   max_sens_spec     0.588         210        210 0.6  
#>  9 9       1         gbm   max_sens_spec     0.462         210        210 0.762
#> 10 10      1         gbm   max_sens_spec     0.512         210        210 0.657
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
#> $data_ens
#> # A tibble: 4,200 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <dbl> <dbl>
#>  1 1      .part1     1         0 0.508
#>  2 6      .part1     1         0 0.406
#>  3 10     .part1     1         0 0.621
#>  4 17     .part1     1         0 0.487
#>  5 18     .part1     1         0 0.380
#>  6 24     .part1     1         0 0.627
#>  7 25     .part1     1         0 0.619
#>  8 28     .part1     1         0 0.550
#>  9 29     .part1     1         0 0.300
#> 10 30     .part1     1         0 0.419
#> # ℹ 4,190 more rows
#> 
# }
```
