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
# \donttest{
require(dplyr)
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

gam_t1 <- fit_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  select_pred = FALSE,
  partition = ".part",
  thr = "max_sens_spec"
)
#> Formula used for model fitting:
#> pr_ab ~ s(aet, k = -1) + s(ppt_jja, k = -1) + s(pH, k = -1) + s(awc, k = -1) + s(depth, k = -1) + landform
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
gam_t1$model
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = -1) + s(ppt_jja, k = -1) + s(pH, k = -1) + 
#>     s(awc, k = -1) + s(depth, k = -1) + landform
#> 
#> Estimated degrees of freedom:
#> 4.06 6.65 7.88 1.00 1.68  total = 36.27 
#> 
#> UBRE score: 0.02901258     
gam_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
gam_t1$performance
#> # A tibble: 1 × 33
#>   model threshold     thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 gam   max_sens_spec     0.530         700        700    0.753  0.114     0.76
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
gam_t1$performance_part
#> # A tibble: 10 × 21
#>    replica partition model threshold     thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>             <dbl>       <int>      <int> <dbl>
#>  1 1       1         gam   max_sens_spec     0.507          70         70 0.771
#>  2 1       2         gam   max_sens_spec     0.702          70         70 0.729
#>  3 1       3         gam   max_sens_spec     0.584          70         70 0.6  
#>  4 1       4         gam   max_sens_spec     0.507          70         70 0.814
#>  5 1       5         gam   max_sens_spec     0.537          70         70 0.729
#>  6 1       6         gam   max_sens_spec     0.383          70         70 0.943
#>  7 1       7         gam   max_sens_spec     0.658          70         70 0.557
#>  8 1       8         gam   max_sens_spec     0.577          70         70 0.757
#>  9 1       9         gam   max_sens_spec     0.555          70         70 0.757
#> 10 1       10        gam   max_sens_spec     0.512          70         70 0.871
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>

# Specifying the formula explicitly
require(mgcv)
#> Loading required package: mgcv
#> Loading required package: nlme
#> 
#> Attaching package: ‘nlme’
#> The following object is masked from ‘package:dplyr’:
#> 
#>     collapse
#> This is mgcv 1.9-4. For overview type '?mgcv'.
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
#> Formula used for model fitting:
#> pr_ab ~ s(aet) + s(ppt_jja) + s(pH) + landform
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

gam_t2$model
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet) + s(ppt_jja) + s(pH) + landform
#> 
#> Estimated degrees of freedom:
#> 4.33 6.70 4.89  total = 30.92 
#> 
#> UBRE score: 0.08202445     
gam_t2$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
gam_t2$performance %>% dplyr::select(ends_with("_mean"))
#> # A tibble: 1 × 14
#>   TPR_mean TNR_mean W_TPR_TNR_mean SORENSEN_mean JACCARD_mean FPB_mean OR_mean
#>      <dbl>    <dbl>          <dbl>         <dbl>        <dbl>    <dbl>   <dbl>
#> 1    0.783    0.703          0.743         0.750        0.602     1.20   0.217
#> # ℹ 7 more variables: TSS_mean <dbl>, KAPPA_mean <dbl>, MCC_mean <dbl>,
#> #   AUC_mean <dbl>, BOYCE_mean <dbl>, CRPS_mean <dbl>, IMAE_mean <dbl>

# Using repeated k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 5)
)
abies2
#> # A tibble: 1,400 × 18
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
#> # ℹ 7 more variables: depth <dbl>, landform <fct>, .part1 <int>, .part2 <int>,
#> #   .part3 <int>, .part4 <int>, .part5 <int>

gam_t3 <- fit_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("ppt_jja", "pH", "awc"),
  predictors_f = c("landform"),
  select_pred = FALSE,
  partition = ".part",
  thr = "max_sens_spec"
)
#> Formula used for model fitting:
#> pr_ab ~ s(ppt_jja, k = -1) + s(pH, k = -1) + s(awc, k = -1) + landform
#> Replica number: 1/5
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
#> Replica number: 2/5
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
#> Replica number: 3/5
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
#> Replica number: 4/5
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
#> Replica number: 5/5
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
gam_t3
#> $model
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = -1) + s(pH, k = -1) + s(awc, k = -1) + 
#>     landform
#> 
#> Estimated degrees of freedom:
#> 6.65 4.04 5.84  total = 31.53 
#> 
#> UBRE score: 0.06214506     
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
#> 1 gam   max_sens_spec     0.554         700        700    0.750 0.0927    0.721
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 25 × 21
#>    replica partition model threshold     thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>             <dbl>       <int>      <int> <dbl>
#>  1 1       1         gam   max_sens_spec     0.653         140        140 0.586
#>  2 1       2         gam   max_sens_spec     0.583         140        140 0.657
#>  3 1       3         gam   max_sens_spec     0.321         140        140 0.957
#>  4 1       4         gam   max_sens_spec     0.595         140        140 0.721
#>  5 1       5         gam   max_sens_spec     0.468         140        140 0.779
#>  6 2       1         gam   max_sens_spec     0.579         140        140 0.743
#>  7 2       2         gam   max_sens_spec     0.554         140        140 0.679
#>  8 2       3         gam   max_sens_spec     0.437         140        140 0.857
#>  9 2       4         gam   max_sens_spec     0.588         140        140 0.664
#> 10 2       5         gam   max_sens_spec     0.535         140        140 0.779
#> # ℹ 15 more rows
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
#> 
#> $data_ens
#> # A tibble: 6,996 × 5
#>    rnames replicates part  pr_ab     pred
#>    <chr>  <chr>      <chr> <dbl>    <dbl>
#>  1 7      .part1     1         0 0.000597
#>  2 10     .part1     1         0 0.511   
#>  3 13     .part1     1         0 0.0585  
#>  4 15     .part1     1         0 0.585   
#>  5 22     .part1     1         0 0.130   
#>  6 33     .part1     1         0 0.0784  
#>  7 35     .part1     1         0 0.493   
#>  8 39     .part1     1         0 0.750   
#>  9 43     .part1     1         0 0.795   
#> 10 44     .part1     1         0 0.705   
#> # ℹ 6,986 more rows
#> 
# }
```
