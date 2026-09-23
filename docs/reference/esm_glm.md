# Fit and validate Generalized Linear Models based on Ensembles of Small Models approach

This function constructs Generalized Linear Models using the Ensembles
of Small Models (ESM) approach (Breiner et al., 2015, 2018).

## Usage

``` r
esm_glm(
  data,
  response,
  predictors,
  partition,
  thr = NULL,
  poly = 0,
  inter_order = 0
)
```

## Arguments

- data:

  data.frame. Database with the response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1).

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). This can only construct models
  with continuous variables and does not allow categorical variables.
  Usage predictors = c("aet", "cwd", "tmin").

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

  - max_jaccard: The threshold at which Jaccard is the highest.

  - max_sorensen: The threshold at which Sorensen is highest.

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

- poly:

  integer \>= 2. If used with values \>= 2 model will use polynomials
  for those continuous variables (i.e. used in predictors argument).
  Default is 0. Because ESM are constructed with few occurrences it is
  recommended no to use polynomials to avoid overfitting.

- inter_order:

  integer \>= 0. The interaction order between explanatory variables.
  Default is 0. Because ESM are constructed with few occurrences it is
  recommended not to use interaction terms.

## Value

A list object with:

- esm_model: A list with "glm" class object from stats package for each
  bivariate model. This object can be used for predicting ensembles of
  small models with
  [`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)
  function.

- predictors: A tibble with variables use for modeling.

- performance: Performance metric (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Those threshold dependent metric are calculated based on the threshold
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
[`esm_gbm`](https://sjevelazco.github.io/flexsdm/reference/esm_gbm.md),
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
  method = c(method = "rep_kfold", folds = 3, replicates = 5)
)
abies2
#> # A tibble: 20 × 18
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
#> # ℹ 7 more variables: depth <dbl>, landform <fct>, .part1 <int>, .part2 <int>,
#> #   .part3 <int>, .part4 <int>, .part5 <int>

# Without threshold specification and with kfold
esm_glm_t1 <- esm_glm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL,
  poly = 0,
  inter_order = 0
)
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   4%
  |                                                                            
  |=====                                                                 |   7%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |============                                                          |  18%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |============================                                          |  39%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |================================                                      |  46%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |========================================                              |  57%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |==========================================================            |  82%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |====================================================================  |  96%
  |                                                                            
  |======================================================================| 100%

esm_glm_t1$esm_model # bivariate model
#> $`0.530555555555555`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet          cwd  
#>    12.61124     -0.01839     -0.01106  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 20.17     AIC: 26.17
#> 
#> $`1`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet         tmin  
#>    -61.7022       0.8884     -98.8127  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 3.623e-09     AIC: 6
#> 
#> $`0.240740740740741`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet      ppt_djf  
#>     3.81661     -0.01077      0.00169  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 23.91     AIC: 29.91
#> 
#> $`0.255555555555556`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet      ppt_jja  
#>     3.71684     -0.00747     -0.04999  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 23.53     AIC: 29.53
#> 
#> $`0.225925925925926`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet           pH  
#>     3.46338     -0.01078      0.15789  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 23.64     AIC: 29.64
#> 
#> $`0.184259259259259`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet          awc  
#>       3.874       -0.009       -4.004  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 23.9  AIC: 29.9
#> 
#> $`0.638888888888889`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          aet        depth  
#>     2.73490     -0.01511      0.02985  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 15.65     AIC: 21.65
#> 
#> $`1`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          cwd         tmin  
#>    255.7081       0.2674    -161.7130  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 1.019e-08     AIC: 6
#> 
#> $`0.389814814814815`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          cwd      ppt_djf  
#>     9.81803     -0.01116     -0.01715  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 22.72     AIC: 28.72
#> 
#> $`0.666666666666667`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          cwd      ppt_jja  
#>    20.03515     -0.02619     -0.53107  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 16.15     AIC: 22.15
#> 
#> $`0.48287037037037`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          cwd        depth  
#>   -4.123825     0.002939     0.024211  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 21.11     AIC: 27.11
#> 
#> $`0.955555555555555`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)         tmin      ppt_djf  
#>     184.029     -194.331        1.404  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 1.158e-08     AIC: 6
#> 
#> $`1`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)         tmin      ppt_jja  
#>      250.78       -82.32        -3.57  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 3.806e-09     AIC: 6
#> 
#> $`0.918981481481481`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)         tmin           pH  
#>      413.24      -109.69       -22.59  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 5.68e-09  AIC: 6
#> 
#> $`0.931944444444444`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)         tmin          awc  
#>      309.09       -83.45     -1100.18  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 4.003e-09     AIC: 6
#> 
#> $`0.881944444444444`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)         tmin        depth  
#>     38.6582     -43.3890       0.6695  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 1.784e-09     AIC: 6
#> 
#> $`0.0925925925925926`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)      ppt_djf      ppt_jja  
#>    1.085941     0.003349    -0.127860  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 24.91     AIC: 30.91
#> 
#> $`0.593518518518519`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)      ppt_djf        depth  
#>    -0.09930     -0.01080      0.02633  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 18.89     AIC: 24.89
#> 
#> $`0.214814814814815`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)      ppt_jja           pH  
#>      1.0179      -0.1114       0.1345  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 24.81     AIC: 30.81
#> 
#> $`0.185185185185185`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)      ppt_jja          awc  
#>     2.20331     -0.08642    -10.81332  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 24.38     AIC: 30.38
#> 
#> $`0.667592592592593`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)      ppt_jja        depth  
#>    -0.21334     -0.18621      0.02613  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 17.32     AIC: 23.32
#> 
#> $`0.0189814814814815`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)           pH          awc  
#>      0.3596       0.3046     -24.1322  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 25.55     AIC: 31.55
#> 
#> $`0.497685185185185`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)           pH        depth  
#>    -0.04891     -0.92171      0.04322  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 16.74     AIC: 22.74
#> 
#> $`0.543518518518519`
#> 
#> Call:  stats::glm(formula = formula1, family = "binomial", data = data)
#> 
#> Coefficients:
#> (Intercept)          awc        depth  
#>    -0.34942    -80.56792      0.06789  
#> 
#> Degrees of Freedom: 19 Total (i.e. Null);  17 Residual
#> Null Deviance:       27.73 
#> Residual Deviance: 11.76     AIC: 17.76
#> 
esm_glm_t1$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_glm_t1$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_glm equal_sens_…     0.472          10         10        1      0        1
#> 2 esm_glm lpt              0.440          10         10        1      0        1
#> 3 esm_glm max_fpb          0.583          10         10        1      0        1
#> 4 esm_glm max_jaccard      0.583          10         10        1      0        1
#> 5 esm_glm max_sens_sp…     0.583          10         10        1      0        1
#> 6 esm_glm max_sorensen     0.583          10         10        1      0        1
#> 7 esm_glm sensitivity      0.614          10         10        1      0        1
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
esm_glm_t1$performancePart
#> NULL
# }
```
