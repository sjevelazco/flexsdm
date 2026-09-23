# Fit and validate Neural Networks models

Fit and validate Neural Networks models

## Usage

``` r
fit_net(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  thr = NULL,
  size = 2,
  decay = 0.1
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
  response, predictors, and predictors_f arguments. Defaul NULL.

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- thr:

  character. Threshold used to get binary suitability values (i.e.
  0,1)., needed for threshold-dependent performance metrics. More than
  one threshold type can be specified. It is necessary to provide a
  vector for this argument. The following threshold criteria are
  available:

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
    specified, the default is 0.9

  If more than one threshold type is used they must be concatenated,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- size:

  numeric. Number of units in the hidden layer. Can be zero if there are
  skip-layer units. Default 2.

- decay:

  numeric. Parameter for weight decay. Default 0.1.

## Value

A list object with:

- model: A "nnet.formula" "nnet" class object from nnet package. This
  object can be used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metric are calculated based on the threshold
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

nnet_t1 <- fit_net(
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

nnet_t1$model
#> a 19-2-1 network with 43 weights
#> inputs: aet ppt_jja pH awc depth landform2 landform3 landform4 landform5 landform6 landform7 landform8 landform9 landform10 landform11 landform12 landform13 landform14 landform15 
#> output(s): pr_ab 
#> options were - entropy fitting  decay=0.1
nnet_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
nnet_t1$performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 net   equal_sens_sp…     0.529         700        700    0.667 0.0670    0.667
#> 2 net   max_sens_spec      0.524         700        700    0.739 0.162     0.66 
#> 3 net   max_sorensen       0.388         700        700    0.886 0.0780    0.467
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
nnet_t1$performance_part
#> # A tibble: 30 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       1         net   max_sorensen       0.388          70         70 0.9  
#>  2 1       1         net   max_sens_spec      0.388          70         70 0.9  
#>  3 1       1         net   equal_sens_sp…     0.512          70         70 0.629
#>  4 1       2         net   max_sorensen       0.192          70         70 0.971
#>  5 1       2         net   max_sens_spec      0.502          70         70 0.743
#>  6 1       2         net   equal_sens_sp…     0.531          70         70 0.7  
#>  7 1       3         net   max_sorensen       0.339          70         70 0.9  
#>  8 1       3         net   max_sens_spec      0.628          70         70 0.614
#>  9 1       3         net   equal_sens_sp…     0.526          70         70 0.7  
#> 10 1       4         net   max_sorensen       0.541          70         70 0.757
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
nnet_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <fct> <dbl>
#>  1 8      .part      1     0     0.131
#>  2 21     .part      1     0     0.283
#>  3 22     .part      1     0     0.151
#>  4 30     .part      1     0     0.501
#>  5 32     .part      1     0     0.279
#>  6 60     .part      1     0     0.519
#>  7 62     .part      1     0     0.425
#>  8 68     .part      1     0     0.646
#>  9 83     .part      1     0     0.282
#> 10 89     .part      1     0     0.479
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

nnet_t2 <- fit_net(
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
nnet_t2
#> $model
#> a 19-2-1 network with 43 weights
#> inputs: aet ppt_jja pH awc depth landform2 landform3 landform4 landform5 landform6 landform7 landform8 landform9 landform10 landform11 landform12 landform13 landform14 landform15 
#> output(s): pr_ab 
#> options were - entropy fitting  decay=0.1
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
#> 1 net   equal_sens_sp…     0.529         700        700    0.666 0.0430    0.666
#> 2 net   max_sens_spec      0.524         700        700    0.736 0.130     0.647
#> 3 net   max_sorensen       0.388         700        700    0.905 0.0469    0.433
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
#>  1 1       1         net   max_sorensen       0.324         210        210 0.9  
#>  2 1       1         net   max_sens_spec      0.539         210        210 0.6  
#>  3 1       1         net   equal_sens_sp…     0.460         210        210 0.652
#>  4 2       1         net   max_sorensen       0.388         210        210 0.886
#>  5 2       1         net   max_sens_spec      0.511         210        210 0.719
#>  6 2       1         net   equal_sens_sp…     0.515         210        210 0.714
#>  7 3       1         net   max_sorensen       0.394         210        210 0.829
#>  8 3       1         net   max_sens_spec      0.420         210        210 0.810
#>  9 3       1         net   equal_sens_sp…     0.499         210        210 0.610
#> 10 4       1         net   max_sorensen       0.344         210        210 0.924
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
#> $data_ens
#> # A tibble: 4,200 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <fct> <dbl>
#>  1 7      .part1     1     0     0.146
#>  2 9      .part1     1     0     0.467
#>  3 13     .part1     1     0     0.440
#>  4 15     .part1     1     0     0.396
#>  5 20     .part1     1     0     0.105
#>  6 23     .part1     1     0     0.305
#>  7 27     .part1     1     0     0.265
#>  8 33     .part1     1     0     0.289
#>  9 39     .part1     1     0     0.290
#> 10 42     .part1     1     0     0.328
#> # ℹ 4,190 more rows
#> 
# }
```
