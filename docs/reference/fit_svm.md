# Fit and validate Support Vector Machine models

Fit and validate Support Vector Machine models

## Usage

``` r
fit_svm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  thr = NULL,
  sigma = "automatic",
  C = 1
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
  response, predictors, and predictors_f arguments

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- thr:

  character. Threshold used to get binary suitability values (i.e. 0,1)
  needed for threshold-dependent performance metrics. More than one
  threshold type can be used. It is necessary to provide a vector for
  this argument. The following threshold criteria are available:

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

- sigma:

  numeric. Inverse kernel width for the Radial Basis kernel function
  "rbfdot". Default "automatic".

- C:

  numeric. Cost of constraints violation, the 'C'-constant of the
  regularization term in the Lagrange formulation. Default 1

## Value

A list object with:

- model: A "ksvm" class object from kernlab package. This object can be
  used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in the argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- data_ens: Predicted suitability for each test partition based on the
  best model. This database is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## Details

This function constructs 'C-svc' classification type and uses Radial
Basis kernel "Gaussian" function (rbfdot). See details details in
[ksvm](https://rdrr.io/pkg/kernlab/man/ksvm.html).

## See also

[`fit_gam`](https://sjevelazco.github.io/flexsdm/reference/fit_gam.md),
[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md),
[`fit_gbm`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md),
[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md),
[`fit_max`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md),
[`fit_net`](https://sjevelazco.github.io/flexsdm/reference/fit_net.md),
and
[`fit_raf`](https://sjevelazco.github.io/flexsdm/reference/fit_raf.md).

## Examples

``` r
data("abies")

# Using k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
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

svm_t1 <- fit_svm(
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
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5

names(svm_t1)
#> [1] "model"            "predictors"       "performance"      "performance_part"
#> [5] "data_ens"        
svm_t1$model
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  0.201598132358874 
#> 
#> Number of Support Vectors : 899 
#> 
#> Objective Function Value : -767.6206 
#> Training error : 0.211429 
#> Probability model included. 
svm_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
svm_t1$performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 svm   equal_sens_sp…     0.541         700        700    0.739 0.0377    0.739
#> 2 svm   max_sens_spec      0.541         700        700    0.8   0.0787    0.717
#> 3 svm   max_sorensen       0.296         700        700    0.88  0.0526    0.613
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
svm_t1$performance_part
#> # A tibble: 15 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       1         svm   max_sorensen       0.511         140        140 0.793
#>  2 1       1         svm   max_sens_spec      0.607         140        140 0.707
#>  3 1       1         svm   equal_sens_sp…     0.557         140        140 0.757
#>  4 1       2         svm   max_sorensen       0.306         140        140 0.921
#>  5 1       2         svm   max_sens_spec      0.598         140        140 0.771
#>  6 1       2         svm   equal_sens_sp…     0.603         140        140 0.757
#>  7 1       3         svm   max_sorensen       0.313         140        140 0.886
#>  8 1       3         svm   max_sens_spec      0.313         140        140 0.886
#>  9 1       3         svm   equal_sens_sp…     0.520         140        140 0.671
#> 10 1       4         svm   max_sorensen       0.309         140        140 0.921
#> 11 1       4         svm   max_sens_spec      0.504         140        140 0.757
#> 12 1       4         svm   equal_sens_sp…     0.495         140        140 0.757
#> 13 1       5         svm   max_sorensen       0.434         140        140 0.879
#> 14 1       5         svm   max_sens_spec      0.434         140        140 0.879
#> 15 1       5         svm   equal_sens_sp…     0.573         140        140 0.75 
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
svm_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab   pred
#>    <chr>  <chr>      <chr> <fct>  <dbl>
#>  1 1      .part      1     0     0.709 
#>  2 5      .part      1     0     0.199 
#>  3 7      .part      1     0     0.0895
#>  4 9      .part      1     0     0.587 
#>  5 17     .part      1     0     0.561 
#>  6 20     .part      1     0     0.214 
#>  7 33     .part      1     0     0.153 
#>  8 41     .part      1     0     0.134 
#>  9 44     .part      1     0     0.478 
#> 10 45     .part      1     0     0.179 
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

svm_t2 <- fit_svm(
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
svm_t2
#> $model
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  0.191996666477365 
#> 
#> Number of Support Vectors : 904 
#> 
#> Objective Function Value : -773.0579 
#> Training error : 0.214286 
#> Probability model included. 
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
#> 1 svm   equal_sens_sp…     0.536         700        700    0.766 0.0156    0.766
#> 2 svm   max_sens_spec      0.544         700        700    0.834 0.0612    0.731
#> 3 svm   max_sorensen       0.286         700        700    0.890 0.0231    0.671
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
#>  1 1       1         svm   max_sorensen       0.417         210        210 0.886
#>  2 1       1         svm   max_sens_spec      0.417         210        210 0.886
#>  3 1       1         svm   equal_sens_sp…     0.576         210        210 0.752
#>  4 2       1         svm   max_sorensen       0.344         210        210 0.938
#>  5 2       1         svm   max_sens_spec      0.555         210        210 0.748
#>  6 2       1         svm   equal_sens_sp…     0.512         210        210 0.776
#>  7 3       1         svm   max_sorensen       0.445         210        210 0.9  
#>  8 3       1         svm   max_sens_spec      0.445         210        210 0.9  
#>  9 3       1         svm   equal_sens_sp…     0.604         210        210 0.752
#> 10 4       1         svm   max_sorensen       0.353         210        210 0.890
#> # ℹ 20 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
#> $data_ens
#> # A tibble: 4,200 × 5
#>    rnames replicates part  pr_ab   pred
#>    <chr>  <chr>      <chr> <fct>  <dbl>
#>  1 6      .part1     1     0     0.544 
#>  2 7      .part1     1     0     0.0617
#>  3 8      .part1     1     0     0.314 
#>  4 10     .part1     1     0     0.320 
#>  5 17     .part1     1     0     0.586 
#>  6 19     .part1     1     0     0.716 
#>  7 20     .part1     1     0     0.175 
#>  8 22     .part1     1     0     0.103 
#>  9 24     .part1     1     0     0.792 
#> 10 27     .part1     1     0     0.133 
#> # ℹ 4,190 more rows
#> 
```
