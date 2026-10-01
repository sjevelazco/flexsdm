# Fit and validate Generalized Additive Models based on Ensembles of Small Models approach

This function constructs Generalized Additive Models using the Ensembles
of Small Models (ESM) approach (Breiner et al., 2015, 2018).

## Usage

``` r
esm_gam(data, response, predictors, partition, thr = NULL, k = 2)
```

## Arguments

- data:

  data.frame. Database with the response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1)

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). This function does not allow
  categorical variables and can only construct models with continuous
  variables. Usage predictors = c("aet", "cwd", "tmin").

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

  If the user wants to include more than one threshold type, it is
  necessary to concatenate threshold types, e.g., thr=c('max_sens_spec',
  'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'),
  or thr=c('max_sens_spec', 'sensitivity'). Function will use all
  thresholds if no threshold is specified

- k:

  integer. The dimension of the basis used to represent the smooth term.
  Default 2. Because ESM was proposed to fit models with little data, we
  recommend using small values of this parameter.

## Value

A list object with:

- esm_model: A list with "gam" class object from mgcv package for each
  bivariate model. This object can be used for predicting an ensemble of
  small models with the
  [`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)
  function.

- predictors: A tibble with variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in the argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

## Details

This method consists of creating bivariate models with all pair-wise
combinations of predictors and perform an ensemble based on the average
of suitability weighted by Somers' D metric (D = 2 x (AUC -0.5)). ESM is
recommended for modeling species with few occurrences. This function
does not allow categorical variables because the use of these types of
variables could be problematic when using with few occurrences. For
further detail see Breiner et al. (2015, 2018).

This function fits GAM using mgvc package, with Binomial distribution
family and thin plate regression spline as a smoothing basis (see
?mgvc::s).

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

[`esm_gau`](https://sjevelazco.github.io/flexsdm/reference/esm_gau.md),
[`esm_gbm`](https://sjevelazco.github.io/flexsdm/reference/esm_gbm.md),
[`esm_glm`](https://sjevelazco.github.io/flexsdm/reference/esm_glm.md),
[`esm_max`](https://sjevelazco.github.io/flexsdm/reference/esm_max.md),
[`esm_net`](https://sjevelazco.github.io/flexsdm/reference/esm_net.md),
and
[`esm_svm`](https://sjevelazco.github.io/flexsdm/reference/esm_svm.md).

## Examples

``` r
# \donttest{
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
  method = c(method = "kfold", folds = 3)
)
abies2
#> # A tibble: 20 × 14
#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc
#>    <int> <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl>  <dbl>
#>  1 12040     0 -308909.  384248.  573.  332.  4.84     521.   48.8   5.63 0.108 
#>  2 10361     0 -254286.  417158.  260.  469.  2.93     151.   15.1   6.20 0.0950
#>  3  9402     0 -286979.  386206.  587.  376.  6.45     333.   15.7   5.5  0.160 
#>  4  9815     0 -291849.  445595.  443.  455.  4.39     332.   19.1   6    0.0700
#>  5 10524     0 -256658.  184438.  355.  568.  5.87     303.   10.6   5.20 0.0800
#>  6  8860     0  121343. -164170.  354.  733.  3.97     182.    9.83  0    0     
#>  7  6431     0  107903. -122968.  461.  578.  4.87     161.    7.66  5.90 0.0900
#>  8 11730     0 -333903.  431238.  561.  364.  6.73     387.   25.2   5.80 0.130 
#>  9   808     0 -150163.  357180.  339.  564.  2.64     220.   15.3   6.40 0.100 
#> 10 11054     0 -293663.  340981.  477.  396.  3.89     332.   26.4   4.60 0.0634
#> 11  2960     1  -49273.  181752.  512.  275.  0.920    319.   17.3   5.92 0.0900
#> 12  3065     1  126907. -198892.  322.  544.  0.700    203.   10.6   5.60 0.110 
#> 13  5527     1  116751. -181089.  261.  537.  0.363    178.    7.43  0    0     
#> 14  4035     1  -31777.  115940.  394.  440.  2.07     298.   11.2   6.01 0.0769
#> 15  4081     1   -5158.   90159.  301.  502.  0.703    203.   14.6   6.11 0.0633
#> 16  3087     1  102151. -143976.  299.  425. -2.08     205.   13.4   3.88 0.110 
#> 17  3495     1  -19586.   89803.  438.  419.  2.13     189.   15.2   6.19 0.0959
#> 18  4441     1   49405.  -60502.  362.  582.  2.42     218.    7.84  5.64 0.0786
#> 19   301     1 -132516.  270845.  367.  196. -2.56     422.   26.3   6.70 0.0300
#> 20  3162     1   59905.  -53634.  319.  626.  1.99     212.    4.50  4.51 0.0396
#> # ℹ 3 more variables: depth <dbl>, landform <fct>, .part <int>

# Without threshold specification and with kfold
esm_gam_t1 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%

esm_gam_t1$esm_model # bivariate model
#> $`0.398148148148148`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(cwd, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.81 1.73  total = 4.54 
#> 
#> UBRE score: 0.2827405     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.0925925925925926`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4818906     
#> 
#> $`0.564814814814815`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.12  total = 3.12 
#> 
#> UBRE score: 0.08148345     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.314814814814815`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4360771     
#> 
#> $`0.875`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.107549     
#> 
#> $`0.553240740740741`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.38  total = 3.38 
#> 
#> UBRE score: 0.3416266     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.888888888888889`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.895833333333333`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.208333333333333`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5455887     
#> 
#> $`0.527777777777778`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.33  total = 3.33 
#> 
#> UBRE score: 0.2354151     
#> 
#> $`0.0092592592592593`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5402802     
#> 
#> $`0.527777777777778`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.36  total = 3.36 
#> 
#> UBRE score: 0.1548222     
#> 
#> $`0.49537037037037`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(pH, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.1369325     
#> 
#> $`0.342592592592593`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(awc, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.1120694     
#> 
esm_gam_t1$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_gam_t1$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_gam equal_sens_…     0.607          10         10        1      0        1
#> 2 esm_gam lpt              0.607          10         10        1      0        1
#> 3 esm_gam max_fpb          0.607          10         10        1      0        1
#> 4 esm_gam max_jaccard      0.607          10         10        1      0        1
#> 5 esm_gam max_sens_sp…     0.607          10         10        1      0        1
#> 6 esm_gam max_sorensen     0.607          10         10        1      0        1
#> 7 esm_gam sensitivity      0.614          10         10        1      0        1
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
esm_gam_t1$performance_part
#> # A tibble: 21 × 21
#>    model replicates part  threshold thr_value n_presences n_absences   TPR   TNR
#>    <chr> <chr>      <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#>  1 esm_… .part      1     max_sore…     0.607           4          4     1     1
#>  2 esm_… .part      1     max_jacc…     0.607           4          4     1     1
#>  3 esm_… .part      1     max_fpb       0.607           4          4     1     1
#>  4 esm_… .part      1     max_sens…     0.607           4          4     1     1
#>  5 esm_… .part      1     equal_se…     0.607           4          4     1     1
#>  6 esm_… .part      1     lpt           0.607           4          4     1     1
#>  7 esm_… .part      1     sensitiv…     0.607           4          4     1     1
#>  8 esm_… .part      2     max_sore…     0.699           3          3     1     1
#>  9 esm_… .part      2     max_jacc…     0.699           3          3     1     1
#> 10 esm_… .part      2     max_fpb       0.699           3          3     1     1
#> # ℹ 11 more rows
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>

# Test with rep_kfold partition
abies2 <- abies2 %>% select(-starts_with("."))

set.seed(10)
abies2 <- part_random(
  data = abies2,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 3, replicates = 2)
)
abies2
#> # A tibble: 20 × 15
#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc
#>    <int> <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl>  <dbl>
#>  1 12040     0 -308909.  384248.  573.  332.  4.84     521.   48.8   5.63 0.108 
#>  2 10361     0 -254286.  417158.  260.  469.  2.93     151.   15.1   6.20 0.0950
#>  3  9402     0 -286979.  386206.  587.  376.  6.45     333.   15.7   5.5  0.160 
#>  4  9815     0 -291849.  445595.  443.  455.  4.39     332.   19.1   6    0.0700
#>  5 10524     0 -256658.  184438.  355.  568.  5.87     303.   10.6   5.20 0.0800
#>  6  8860     0  121343. -164170.  354.  733.  3.97     182.    9.83  0    0     
#>  7  6431     0  107903. -122968.  461.  578.  4.87     161.    7.66  5.90 0.0900
#>  8 11730     0 -333903.  431238.  561.  364.  6.73     387.   25.2   5.80 0.130 
#>  9   808     0 -150163.  357180.  339.  564.  2.64     220.   15.3   6.40 0.100 
#> 10 11054     0 -293663.  340981.  477.  396.  3.89     332.   26.4   4.60 0.0634
#> 11  2960     1  -49273.  181752.  512.  275.  0.920    319.   17.3   5.92 0.0900
#> 12  3065     1  126907. -198892.  322.  544.  0.700    203.   10.6   5.60 0.110 
#> 13  5527     1  116751. -181089.  261.  537.  0.363    178.    7.43  0    0     
#> 14  4035     1  -31777.  115940.  394.  440.  2.07     298.   11.2   6.01 0.0769
#> 15  4081     1   -5158.   90159.  301.  502.  0.703    203.   14.6   6.11 0.0633
#> 16  3087     1  102151. -143976.  299.  425. -2.08     205.   13.4   3.88 0.110 
#> 17  3495     1  -19586.   89803.  438.  419.  2.13     189.   15.2   6.19 0.0959
#> 18  4441     1   49405.  -60502.  362.  582.  2.42     218.    7.84  5.64 0.0786
#> 19   301     1 -132516.  270845.  367.  196. -2.56     422.   26.3   6.70 0.0300
#> 20  3162     1   59905.  -53634.  319.  626.  1.99     212.    4.50  4.51 0.0396
#> # ℹ 4 more variables: depth <dbl>, landform <fct>, .part1 <int>, .part2 <int>

esm_gam_t2 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
esm_gam_t2$esm_model # bivariate model
#> $`0.236111111111111`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(cwd, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.81 1.73  total = 4.54 
#> 
#> UBRE score: 0.2827405     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.268518518518519`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4765413     
#> 
#> $`0.0833333333333333`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4818906     
#> 
#> $`0.273148148148148`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4948666     
#> 
#> $`0.546296296296296`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.12  total = 3.12 
#> 
#> UBRE score: 0.08148345     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.342592592592593`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4360771     
#> 
#> $`0.664351851851852`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.107549     
#> 
#> $`0.141203703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5582199     
#> 
#> $`0.392361111111111`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.38  total = 3.38 
#> 
#> UBRE score: 0.3416266     
#> 
#> $`0.944444444444444`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.833333333333333`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.944444444444444`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.888888888888889`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.115740740740741`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5455887     
#> 
#> $`0.0949074074074074`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.619206     
#> 
#> $`0.384259259259259`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.33  total = 3.33 
#> 
#> UBRE score: 0.2354151     
#> 
#> $`0.226851851851852`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5402802     
#> 
#> $`0.115740740740741`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5189453     
#> 
#> $`0.465277777777778`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.36  total = 3.36 
#> 
#> UBRE score: 0.1548222     
#> 
#> $`0.527777777777778`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(pH, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.1369325     
#> 
#> $`0.493055555555556`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(awc, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.1120694     
#> 
esm_gam_t2$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_gam_t2$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_gam equal_sens_…     0.656          10         10    0.958  0.102    0.958
#> 2 esm_gam lpt              0.584          10         10    1      0        0.958
#> 3 esm_gam max_fpb          0.584          10         10    1      0        0.958
#> 4 esm_gam max_jaccard      0.584          10         10    1      0        0.958
#> 5 esm_gam max_sens_sp…     0.656          10         10    0.958  0.102    1    
#> 6 esm_gam max_sorensen     0.584          10         10    1      0        0.958
#> 7 esm_gam sensitivity      0.680          10         10    1      0        0.958
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>

# Test with other bootstrap
abies2 <- abies2 %>% select(-starts_with("."))

set.seed(10)
abies2 <- part_random(
  data = abies2,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 3, proportion = 0.7)
)
abies2
#> # A tibble: 20 × 16
#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc
#>    <int> <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl>  <dbl>
#>  1 12040     0 -308909.  384248.  573.  332.  4.84     521.   48.8   5.63 0.108 
#>  2 10361     0 -254286.  417158.  260.  469.  2.93     151.   15.1   6.20 0.0950
#>  3  9402     0 -286979.  386206.  587.  376.  6.45     333.   15.7   5.5  0.160 
#>  4  9815     0 -291849.  445595.  443.  455.  4.39     332.   19.1   6    0.0700
#>  5 10524     0 -256658.  184438.  355.  568.  5.87     303.   10.6   5.20 0.0800
#>  6  8860     0  121343. -164170.  354.  733.  3.97     182.    9.83  0    0     
#>  7  6431     0  107903. -122968.  461.  578.  4.87     161.    7.66  5.90 0.0900
#>  8 11730     0 -333903.  431238.  561.  364.  6.73     387.   25.2   5.80 0.130 
#>  9   808     0 -150163.  357180.  339.  564.  2.64     220.   15.3   6.40 0.100 
#> 10 11054     0 -293663.  340981.  477.  396.  3.89     332.   26.4   4.60 0.0634
#> 11  2960     1  -49273.  181752.  512.  275.  0.920    319.   17.3   5.92 0.0900
#> 12  3065     1  126907. -198892.  322.  544.  0.700    203.   10.6   5.60 0.110 
#> 13  5527     1  116751. -181089.  261.  537.  0.363    178.    7.43  0    0     
#> 14  4035     1  -31777.  115940.  394.  440.  2.07     298.   11.2   6.01 0.0769
#> 15  4081     1   -5158.   90159.  301.  502.  0.703    203.   14.6   6.11 0.0633
#> 16  3087     1  102151. -143976.  299.  425. -2.08     205.   13.4   3.88 0.110 
#> 17  3495     1  -19586.   89803.  438.  419.  2.13     189.   15.2   6.19 0.0959
#> 18  4441     1   49405.  -60502.  362.  582.  2.42     218.    7.84  5.64 0.0786
#> 19   301     1 -132516.  270845.  367.  196. -2.56     422.   26.3   6.70 0.0300
#> 20  3162     1   59905.  -53634.  319.  626.  1.99     212.    4.50  4.51 0.0396
#> # ℹ 5 more variables: depth <dbl>, landform <fct>, .part1 <chr>, .part2 <chr>,
#> #   .part3 <chr>

esm_gam_t3 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
esm_gam_t3$esm_model # bivariate model
#> $`0.555555555555556`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(cwd, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.81 1.73  total = 4.54 
#> 
#> UBRE score: 0.2827405     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4954043     
#> 
#> $`0.925925925925926`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4765413     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4818906     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4948666     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(aet, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.12  total = 3.12 
#> 
#> UBRE score: 0.08148345     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(tmin, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.185185185185185`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.4360771     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.107549     
#> 
#> $`0.555555555555556`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5582199     
#> 
#> $`0.925925925925926`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(cwd, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.38  total = 3.38 
#> 
#> UBRE score: 0.3416266     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_djf, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(tmin, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.7     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(ppt_jja, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5455887     
#> 
#> $`0.037037037037037`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.619206     
#> 
#> $`0.481481481481481`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5759198     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_djf, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.33  total = 3.33 
#> 
#> UBRE score: 0.2354151     
#> 
#> $`0.703703703703704`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(pH, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5402802     
#> 
#> $`0.777777777777778`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5189453     
#> 
#> $`0.925925925925926`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(ppt_jja, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1.00 1.36  total = 3.36 
#> 
#> UBRE score: 0.1548222     
#> 
#> $`0.555555555555556`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(pH, k = 2) + s(awc, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.5776162     
#> 
#> $`0.925925925925926`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(pH, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.1369325     
#> 
#> $`1`
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pr_ab ~ s(awc, k = 2) + s(depth, k = 2)
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: -0.1120694     
#> 
esm_gam_t3$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_gam_t3$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_gam equal_sens_…     0.686          10         10        1      0        1
#> 2 esm_gam lpt              0.686          10         10        1      0        1
#> 3 esm_gam max_fpb          0.686          10         10        1      0        1
#> 4 esm_gam max_jaccard      0.686          10         10        1      0        1
#> 5 esm_gam max_sens_sp…     0.686          10         10        1      0        1
#> 6 esm_gam max_sorensen     0.686          10         10        1      0        1
#> 7 esm_gam sensitivity      0.686          10         10        1      0        1
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
# }
```
