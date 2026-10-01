# Fit and validate Support Vector Machine models based on Ensembles of Small of Models approach

This function constructs Support Vector Machine models using the
Ensembles of Small Models (ESM) approach (Breiner et al., 2015, 2018).

## Usage

``` r
esm_svm(
  data,
  response,
  predictors,
  partition,
  thr = NULL,
  sigma = "automatic",
  C = 1
)
```

## Arguments

- data:

  data.frame. Database with the response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1).

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). This function only can
  construct models with continuous variables and does not allow
  categorical variables. Usage predictors = c("aet", "cwd", "tmin").

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
    specified, the default is 0.9

  If the user wants to include more than one threshold type, it is
  necessary concatenate threshold types, e.g., thr=c('max_sens_spec',
  'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'),
  or thr=c('max_sens_spec', 'sensitivity'). Function will use all
  thresholds if no threshold is specified

- sigma:

  numeric. Inverse kernel width for the Radial Basis kernel function
  "rbfdot". Default "automatic".

- C:

  numeric. Cost of constraints violation, the 'C' constant of the
  regularization term in the Lagrange formulation. Default 1

## Value

A list object with:

- esm_model: A list with "ksvm" class object from ksvm package for each
  bivariate model. This object can be used for predicting ensemble of
  small model with
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
variables could be problematic when using with few occurrences. Further
detail see Breiner et al. (2015, 2018). This function constructs 'C-svc'
classification type and uses Radial Basis kernel "Gaussian" function
(rbfdot). See details in
[ksvm](https://rdrr.io/pkg/kernlab/man/ksvm.html)

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
[`esm_glm`](https://sjevelazco.github.io/flexsdm/reference/esm_glm.md),
[`esm_max`](https://sjevelazco.github.io/flexsdm/reference/esm_max.md),,
and
[`esm_net`](https://sjevelazco.github.io/flexsdm/reference/esm_net.md).

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
esm_svm_t1 <- esm_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c(
    "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
    "pH", "awc", "depth"
  ),
  partition = ".part",
  thr = NULL
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%

esm_svm_t1$esm_model # bivariate model
#> $`0.867592592592593`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  3.21964131577722 
#> 
#> Number of Support Vectors : 19 
#> 
#> Objective Function Value : -7.6598 
#> Training error : 0.05 
#> Probability model included. 
#> 
#> $`0.74537037037037`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  2.36704511225715 
#> 
#> Number of Support Vectors : 17 
#> 
#> Objective Function Value : -7.5493 
#> Training error : 0.05 
#> Probability model included. 
#> 
#> $`0.863888888888889`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  1.10350006711814 
#> 
#> Number of Support Vectors : 16 
#> 
#> Objective Function Value : -7.8686 
#> Training error : 0.05 
#> Probability model included. 
#> 
#> $`0.219907407407407`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  3.15364895100329 
#> 
#> Number of Support Vectors : 19 
#> 
#> Objective Function Value : -14.3366 
#> Training error : 0.15 
#> Probability model included. 
#> 
#> $`0.303703703703704`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  1.9059427923859 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -14.2805 
#> Training error : 0.2 
#> Probability model included. 
#> 
#> $`0.436111111111111`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  0.890015879420413 
#> 
#> Number of Support Vectors : 16 
#> 
#> Objective Function Value : -8.9456 
#> Training error : 0.1 
#> Probability model included. 
#> 
#> $`0.928703703703704`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  1.21113352969592 
#> 
#> Number of Support Vectors : 15 
#> 
#> Objective Function Value : -7.2926 
#> Training error : 0.05 
#> Probability model included. 
#> 
#> $`0.968518518518519`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  6.69160658225889 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -8.2629 
#> Training error : 0 
#> Probability model included. 
#> 
#> $`0.878240740740741`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  24.6256479156547 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -7.7436 
#> Training error : 0 
#> Probability model included. 
#> 
#> $`0.874074074074074`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  2.76404088910699 
#> 
#> Number of Support Vectors : 17 
#> 
#> Objective Function Value : -8.1925 
#> Training error : 0 
#> Probability model included. 
#> 
#> $`1`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  10.1022680516999 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -7.4825 
#> Training error : 0 
#> Probability model included. 
#> 
#> $`0.169444444444445`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  1.14332635493013 
#> 
#> Number of Support Vectors : 16 
#> 
#> Objective Function Value : -9.7946 
#> Training error : 0.15 
#> Probability model included. 
#> 
#> $`0.239351851851852`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  6.48807292162547 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -15.1365 
#> Training error : 0.15 
#> Probability model included. 
#> 
#> $`0.0092592592592593`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  32.0153706984081 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -12.1117 
#> Training error : 0.05 
#> Probability model included. 
#> 
#> $`0.381481481481482`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  0.495342417602049 
#> 
#> Number of Support Vectors : 15 
#> 
#> Objective Function Value : -10.3477 
#> Training error : 0.15 
#> Probability model included. 
#> 
#> $`0.24212962962963`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  17.2551055210103 
#> 
#> Number of Support Vectors : 20 
#> 
#> Objective Function Value : -15.1682 
#> Training error : 0.2 
#> Probability model included. 
#> 
#> $`0.379166666666667`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  9.36572774490608 
#> 
#> Number of Support Vectors : 19 
#> 
#> Objective Function Value : -9.8635 
#> Training error : 0.1 
#> Probability model included. 
#> 
#> $`0.529166666666667`
#> Support Vector Machine object of class "ksvm" 
#> 
#> SV type: C-svc  (classification) 
#>  parameter : cost C = 1 
#> 
#> Gaussian Radial Basis kernel function. 
#>  Hyperparameter : sigma =  3.23120038498244 
#> 
#> Number of Support Vectors : 19 
#> 
#> Objective Function Value : -8.8474 
#> Training error : 0.05 
#> Probability model included. 
#> 
esm_svm_t1$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_svm_t1$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_svm equal_sens_…     0.564          10         10    0.978 0.0861    0.978
#> 2 esm_svm lpt              0.504          10         10    1     0         0.978
#> 3 esm_svm max_fpb          0.557          10         10    1     0         0.978
#> 4 esm_svm max_jaccard      0.557          10         10    1     0         0.978
#> 5 esm_svm max_sens_sp…     0.565          10         10    0.978 0.0861    1    
#> 6 esm_svm max_sorensen     0.557          10         10    1     0         0.978
#> 7 esm_svm sensitivity      0.588          10         10    1     0         0.978
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
esm_svm_t1$performance_part
#> # A tibble: 105 × 21
#>    model replicates part  threshold thr_value n_presences n_absences   TPR   TNR
#>    <chr> <chr>      <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#>  1 esm_… .part1     1     max_sore…     0.573           4          4     1     1
#>  2 esm_… .part1     1     max_jacc…     0.573           4          4     1     1
#>  3 esm_… .part1     1     max_fpb       0.573           4          4     1     1
#>  4 esm_… .part1     1     max_sens…     0.573           4          4     1     1
#>  5 esm_… .part1     1     equal_se…     0.573           4          4     1     1
#>  6 esm_… .part1     1     lpt           0.573           4          4     1     1
#>  7 esm_… .part1     1     sensitiv…     0.573           4          4     1     1
#>  8 esm_… .part1     2     max_sore…     0.712           3          3     1     1
#>  9 esm_… .part1     2     max_jacc…     0.712           3          3     1     1
#> 10 esm_… .part1     2     max_fpb       0.712           3          3     1     1
#> # ℹ 95 more rows
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
# }
```
