# Fit and validate Maximum Entropy Models based on Ensemble of Small of Model approach

This function constructs Maxent Models using the Ensemble of Small Model
(ESM) approach (Breiner et al., 2015, 2018).

## Usage

``` r
esm_max(
  data,
  response,
  predictors,
  partition,
  thr = NULL,
  background = NULL,
  clamp = TRUE,
  classes = "default",
  pred_type = "cloglog",
  regmult = 2.5
)
```

## Arguments

- data:

  data.frame. Database with the response (0,1) and predictors values.

- response:

  character. Column name with species absence-presence data (0,1)

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). This function can only
  construct models with continuous variables, and does not allow
  categorical variables Usage predictors = c("aet", "cwd", "tmin").

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
    refers to sensitivity value. If no sensitivity value is specified,
    the default is 0.9

  If the user wants to include more than one threshold type, it is
  necessary concatenate threshold types, e.g., thr=c('max_sens_spec',
  'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'),
  or thr=c('max_sens_spec', 'sensitivity'). Function will use all
  thresholds if no threshold is specified.

- background:

  data.frame. Database with response column only with 0 and predictors
  variables. All column names must be consistent with data. Default NULL

- clamp:

  logical. It is set with TRUE, predictors and features are restricted
  to the range seen during model training.

- classes:

  character. A single feature of any combinations of them. Features are
  symbolized by letters: l (linear), q (quadratic), h (hinge), p
  (product), and t (threshold). Usage classes = "lpq". Default "default"
  (see details).

- pred_type:

  character. Type of response required available "link", "exponential",
  "cloglog" and "logistic". Default "cloglog"

- regmult:

  numeric. A constant to adjust regularization. Because ESM are used for
  modeling species with few records default value is 2.5

## Value

A list object with:

- esm_model: A list with "maxnet" class object from maxnet package for
  each bivariate model. This object can be used for predicting ensembles
  of small models with
  [`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)
  function.

- predictors: A tibble with variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Those threshold dependent metric are calculated based on the threshold
  specified in the argument.

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
further detail see Breiner et al. (2015, 2018). This function use a
default regularization multiplier equal to 2.5 (see Breiner et al.,
2018)

When the argument “classes” is set as default MaxEnt will use different
features combination depending of the number of presences (np) with the
follow rule: if np \< 10 classes = "l", if np between 10 and 15 classes
= "lq", if np between 15 and 80 classes = "lqh", and if np \>= 80
classes = "lqph"

When presence-absence (or presence-pseudo-absence) data are used in data
argument in addition to background points, the function will fit models
with presences and background points and validate with presences and
absences. This procedure makes maxent comparable to other
presences-absences models (e.g., random forest, support vector machine).
If only presences and background points data are used, function will fit
and validate model with presences and background data. If only
presence-absences are used in data argument and without background,
function will fit model with the specified data (not recommended).

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
[`esm_net`](https://sjevelazco.github.io/flexsdm/reference/esm_net.md),
and
[`esm_svm`](https://sjevelazco.github.io/flexsdm/reference/esm_svm.md).

## Examples

``` r
# \donttest{
data("abies")
data("backg")
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

set.seed(10)
backg2 <- backg %>%
  na.omit() %>%
  group_by(pr_ab) %>%
  dplyr::slice_sample(n = 100) %>%
  group_by()

backg2 <- part_random(
  data = backg2,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 3, replicates = 2)
)
backg2
#> # A tibble: 100 × 15
#>    pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja     pH      awc
#>    <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl>  <dbl>    <dbl>
#>  1     0  -23361.  129722.  448.  257.  0.683   285.   19.6   0.0230 0.000356
#>  2     0   54669. -311728.  149. 1328. 10.9      24.7   0.710 8.20   0.0900  
#>  3     0 -162141.  -59278.  332.  863.  9.56    110.    1.92  6.34   0.138   
#>  4     0  185349. -426208.  380. 1105. 11.6     105.    3.04  0.588  0.00995 
#>  5     0  125949. -180238.  220.  303. -4.75    216.    9.95  0      0       
#>  6     0  -27411. -378148.  322. 1047.  8.13     87.8   0.717 7.73   0.141   
#>  7     0   17409.  -99508.  340. 1003.  8.15     93.7   1.49  6.5    0.0708  
#>  8     0 -216411.  432662.  259.  657.  4.44    110.   16.3   5.37   0.0704  
#>  9     0  311709. -456718.  271.  877.  6.34     83.0   8.45  0      0       
#> 10     0 -211821.  419972.  369.  633.  3.36     73.2  16.6   6.12   0.158   
#> # ℹ 90 more rows
#> # ℹ 5 more variables: depth <dbl>, percent_clay <dbl>, landform <fct>,
#> #   .part1 <int>, .part2 <int>

# Without threshold specification and with kfold
esm_max_t1 <- esm_max(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL,
  background = backg2,
  clamp = TRUE,
  classes = "default",
  pred_type = "cloglog",
  regmult = 1
)
#>   |                                                                              |                                                                      |   0%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==                                                                    |   4%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |=====                                                                 |   7%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |========                                                              |  11%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==========                                                            |  14%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |============                                                          |  18%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |===============                                                       |  21%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==================                                                    |  25%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |====================                                                  |  29%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |======================                                                |  32%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |=========================                                             |  36%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |============================                                          |  39%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==============================                                        |  43%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |================================                                      |  46%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |===================================                                   |  50%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |======================================                                |  54%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |========================================                              |  57%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==========================================                            |  61%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |=============================================                         |  64%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |================================================                      |  68%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==================================================                    |  71%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |====================================================                  |  75%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |=======================================================               |  79%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==========================================================            |  82%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |============================================================          |  86%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |==============================================================        |  89%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |=================================================================     |  93%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |====================================================================  |  96%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |======================================================================| 100%

esm_max_t1$esm_model # bivariate model
#> $`0.141203703703704`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  97880
#> 2    0  0.00  93460
#> 3    0  0.00  89230
#> 4    0  0.00  85190
#> 5    0  0.00  81340
#> 6    0  0.00  77660
#> 7    0  0.00  74150
#> 8    0  0.00  70790
#> 9    0  0.00  67590
#> 10   0  0.00  64540
#> 11   0  0.00  61620
#> 12   0  0.00  58830
#> 13   0  0.00  56170
#> 14   0  0.00  53630
#> 15   0  0.00  51200
#> 16   0  0.00  48890
#> 17   0  0.00  46680
#> 18   0  0.00  44570
#> 19   0  0.00  42550
#> 20   0  0.00  40630
#> 21   0  0.00  38790
#> 22   0  0.00  37030
#> 23   0  0.00  35360
#> 24   0  0.00  33760
#> 25   0  0.00  32230
#> 26   0  0.00  30770
#> 27   0  0.00  29380
#> 28   0  0.00  28050
#> 29   0  0.00  26780
#> 30   0  0.00  25570
#> 31   0  0.00  24420
#> 32   0  0.00  23310
#> 33   0  0.00  22260
#> 34   0  0.00  21250
#> 35   0  0.00  20290
#> 36   0  0.00  19370
#> 37   0  0.00  18500
#> 38   0  0.00  17660
#> 39   0  0.00  16860
#> 40   0  0.00  16100
#> 41   0  0.00  15370
#> 42   0  0.00  14680
#> 43   0  0.00  14010
#> 44   0  0.00  13380
#> 45   0  0.00  12770
#> 46   0  0.00  12200
#> 47   0  0.00  11640
#> 48   0  0.00  11120
#> 49   0  0.00  10610
#> 50   0  0.00  10130
#> 51   0  0.00   9676
#> 52   0  0.00   9238
#> 53   0  0.00   8820
#> 54   0  0.00   8421
#> 55   0  0.00   8040
#> 56   0  0.00   7677
#> 57   0  0.00   7330
#> 58   0  0.00   6998
#> 59   0  0.00   6682
#> 60   0  0.00   6379
#> 61   0  0.00   6091
#> 62   0  0.00   5815
#> 63   0  0.00   5552
#> 64   0  0.00   5301
#> 65   0  0.00   5061
#> 66   0  0.00   4833
#> 67   0  0.00   4614
#> 68   0  0.00   4405
#> 69   0  0.00   4206
#> 70   0  0.00   4016
#> 71   0  0.00   3834
#> 72   0  0.00   3661
#> 73   0  0.00   3495
#> 74   0  0.00   3337
#> 75   0  0.00   3186
#> 76   0  0.00   3042
#> 77   0  0.00   2904
#> 78   0  0.00   2773
#> 79   0  0.00   2648
#> 80   0  0.00   2528
#> 81   0  0.00   2414
#> 82   0  0.00   2304
#> 83   0  0.00   2200
#> 84   0  0.00   2101
#> 85   0  0.00   2006
#> 86   0  0.00   1915
#> 87   0  0.00   1828
#> 88   0  0.00   1746
#> 89   0  0.00   1667
#> 90   0  0.00   1591
#> 91   0  0.00   1519
#> 92   0  0.00   1451
#> 93   0  0.00   1385
#> 94   0  0.00   1322
#> 95   0  0.00   1263
#> 96   0  0.00   1205
#> 97   0  0.00   1151
#> 98   0  0.00   1099
#> 99   0  0.00   1049
#> 100  0  0.00   1002
#> 101  0  0.00    956
#> 102  0  0.00    913
#> 103  0  0.00    872
#> 104  0  0.00    832
#> 105  0  0.00    795
#> 106  0  0.00    759
#> 107  0  0.00    724
#> 108  0  0.00    692
#> 109  0  0.00    660
#> 110  0  0.00    631
#> 111  0  0.00    602
#> 112  0  0.00    575
#> 113  0  0.00    549
#> 114  0  0.00    524
#> 115  0  0.00    500
#> 116  0  0.00    478
#> 117  0  0.00    456
#> 118  0  0.00    436
#> 119  0  0.00    416
#> 120  0  0.00    397
#> 121  0  0.00    379
#> 122  0  0.00    362
#> 123  0  0.00    346
#> 124  0  0.00    330
#> 125  0  0.00    315
#> 126  0  0.00    301
#> 127  0  0.00    287
#> 128  0  0.00    274
#> 129  0  0.00    262
#> 130  0  0.00    250
#> 131  0  0.00    239
#> 132  0  0.00    228
#> 133  0  0.00    218
#> 134  0  0.00    208
#> 135  0  0.00    198
#> 136  0  0.00    189
#> 137  0  0.00    181
#> 138  0  0.00    173
#> 139  0  0.00    165
#> 140  1  0.10    157
#> 141  1  0.77    150
#> 142  1  1.39    143
#> 143  1  1.97    137
#> 144  1  2.52    131
#> 145  1  3.03    125
#> 146  1  3.51    119
#> 147  1  3.97    114
#> 148  1  4.39    109
#> 149  1  4.80    104
#> 150  1  5.18     99
#> 151  1  5.54     95
#> 152  1  5.87     90
#> 153  1  6.20     86
#> 154  1  6.50     82
#> 155  1  6.79     79
#> 156  1  7.06     75
#> 157  1  7.32     72
#> 158  1  7.56     68
#> 159  1  7.79     65
#> 160  1  8.01     62
#> 161  1  8.21     60
#> 162  1  8.41     57
#> 163  1  8.59     54
#> 164  1  8.77     52
#> 165  1  8.93     49
#> 166  1  9.09     47
#> 167  1  9.24     45
#> 168  1  9.38     43
#> 169  1  9.51     41
#> 170  1  9.64     39
#> 171  1  9.75     37
#> 172  1  9.86     36
#> 173  1  9.97     34
#> 174  1 10.06     33
#> 175  1 10.16     31
#> 176  1 10.24     30
#> 177  1 10.33     28
#> 178  1 10.40     27
#> 179  1 10.47     26
#> 180  1 10.54     25
#> 181  1 10.61     24
#> 182  1 10.66     23
#> 183  1 10.72     22
#> 184  1 10.77     21
#> 185  1 10.82     20
#> 186  1 10.86     19
#> 187  1 10.91     18
#> 188  1 10.95     17
#> 189  1 10.98     16
#> 190  1 11.02     16
#> 191  1 11.05     15
#> 192  2 11.11     14
#> 193  2 11.18     14
#> 194  2 11.24     13
#> 195  2 11.29     12
#> 196  2 11.34     12
#> 197  2 11.39     11
#> 198  2 11.43     11
#> 199  2 11.47     10
#> 200  2 11.51     10
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  33480
#> 2    0  0.00  31970
#> 3    0  0.00  30520
#> 4    0  0.00  29140
#> 5    0  0.00  27820
#> 6    0  0.00  26570
#> 7    0  0.00  25360
#> 8    0  0.00  24220
#> 9    0  0.00  23120
#> 10   0  0.00  22080
#> 11   0  0.00  21080
#> 12   0  0.00  20120
#> 13   0  0.00  19210
#> 14   0  0.00  18350
#> 15   0  0.00  17520
#> 16   0  0.00  16720
#> 17   0  0.00  15970
#> 18   0  0.00  15240
#> 19   0  0.00  14560
#> 20   0  0.00  13900
#> 21   0  0.00  13270
#> 22   0  0.00  12670
#> 23   0  0.00  12100
#> 24   0  0.00  11550
#> 25   0  0.00  11030
#> 26   0  0.00  10530
#> 27   0  0.00  10050
#> 28   0  0.00   9597
#> 29   0  0.00   9163
#> 30   0  0.00   8748
#> 31   0  0.00   8353
#> 32   0  0.00   7975
#> 33   0  0.00   7614
#> 34   0  0.00   7270
#> 35   0  0.00   6941
#> 36   0  0.00   6627
#> 37   0  0.00   6327
#> 38   0  0.00   6041
#> 39   0  0.00   5768
#> 40   0  0.00   5507
#> 41   0  0.00   5258
#> 42   0  0.00   5020
#> 43   0  0.00   4793
#> 44   0  0.00   4576
#> 45   0  0.00   4369
#> 46   0  0.00   4172
#> 47   0  0.00   3983
#> 48   0  0.00   3803
#> 49   0  0.00   3631
#> 50   0  0.00   3467
#> 51   0  0.00   3310
#> 52   0  0.00   3160
#> 53   0  0.00   3017
#> 54   0  0.00   2881
#> 55   0  0.00   2750
#> 56   0  0.00   2626
#> 57   0  0.00   2507
#> 58   0  0.00   2394
#> 59   0  0.00   2286
#> 60   0  0.00   2182
#> 61   0  0.00   2084
#> 62   0  0.00   1989
#> 63   0  0.00   1899
#> 64   0  0.00   1813
#> 65   0  0.00   1731
#> 66   0  0.00   1653
#> 67   0  0.00   1578
#> 68   0  0.00   1507
#> 69   0  0.00   1439
#> 70   0  0.00   1374
#> 71   0  0.00   1312
#> 72   0  0.00   1252
#> 73   0  0.00   1196
#> 74   0  0.00   1142
#> 75   0  0.00   1090
#> 76   0  0.00   1041
#> 77   0  0.00    994
#> 78   0  0.00    949
#> 79   0  0.00    906
#> 80   0  0.00    865
#> 81   0  0.00    826
#> 82   0  0.00    788
#> 83   0  0.00    753
#> 84   0  0.00    719
#> 85   0  0.00    686
#> 86   0  0.00    655
#> 87   0  0.00    625
#> 88   0  0.00    597
#> 89   0  0.00    570
#> 90   0  0.00    544
#> 91   0  0.00    520
#> 92   0  0.00    496
#> 93   0  0.00    474
#> 94   0  0.00    452
#> 95   0  0.00    432
#> 96   0  0.00    412
#> 97   0  0.00    394
#> 98   0  0.00    376
#> 99   0  0.00    359
#> 100  0  0.00    343
#> 101  0  0.00    327
#> 102  0  0.00    312
#> 103  0  0.00    298
#> 104  0  0.00    285
#> 105  1  0.55    272
#> 106  1  1.35    260
#> 107  1  2.09    248
#> 108  1  2.78    237
#> 109  1  3.42    226
#> 110  1  4.02    216
#> 111  1  4.59    206
#> 112  1  5.12    197
#> 113  1  5.62    188
#> 114  1  6.10    179
#> 115  1  6.54    171
#> 116  1  6.96    163
#> 117  1  7.36    156
#> 118  1  7.73    149
#> 119  1  8.09    142
#> 120  1  8.43    136
#> 121  1  8.76    130
#> 122  1  9.07    124
#> 123  1  9.36    118
#> 124  1  9.64    113
#> 125  1  9.91    108
#> 126  1 10.17    103
#> 127  1 10.41     98
#> 128  1 10.65     94
#> 129  1 10.87     90
#> 130  1 11.08     85
#> 131  1 11.29     82
#> 132  1 11.49     78
#> 133  1 11.68     74
#> 134  1 11.86     71
#> 135  1 12.03     68
#> 136  1 12.20     65
#> 137  1 12.36     62
#> 138  1 12.51     59
#> 139  1 12.66     56
#> 140  1 12.80     54
#> 141  1 12.94     51
#> 142  1 13.07     49
#> 143  1 13.20     47
#> 144  1 13.32     45
#> 145  1 13.44     43
#> 146  1 13.55     41
#> 147  1 13.66     39
#> 148  1 13.76     37
#> 149  1 13.86     35
#> 150  1 13.96     34
#> 151  1 14.05     32
#> 152  1 14.14     31
#> 153  1 14.22     29
#> 154  1 14.31     28
#> 155  1 14.38     27
#> 156  1 14.46     26
#> 157  1 14.53     24
#> 158  1 14.60     23
#> 159  1 14.66     22
#> 160  1 14.73     21
#> 161  1 14.79     20
#> 162  1 14.85     19
#> 163  1 14.90     19
#> 164  1 14.95     18
#> 165  1 15.00     17
#> 166  1 15.05     16
#> 167  1 15.09     15
#> 168  1 15.14     15
#> 169  1 15.18     14
#> 170  1 15.22     13
#> 171  1 15.25     13
#> 172  1 15.29     12
#> 173  1 15.32     12
#> 174  1 15.35     11
#> 175  1 15.38     11
#> 176  1 15.41     10
#> 177  1 15.43     10
#> 178  1 15.46      9
#> 179  1 15.48      9
#> 180  1 15.50      8
#> 181  2 15.59      8
#> 182  2 15.82      8
#> 183  2 16.02      7
#> 184  2 16.21      7
#> 185  2 16.38      7
#> 186  2 16.53      6
#> 187  2 16.66      6
#> 188  2 16.79      6
#> 189  2 16.90      6
#> 190  2 17.00      5
#> 191  2 17.09      5
#> 192  2 17.17      5
#> 193  2 17.25      5
#> 194  2 17.31      4
#> 195  2 17.37      4
#> 196  2 17.43      4
#> 197  2 17.48      4
#> 198  2 17.53      4
#> 199  2 17.57      4
#> 200  2 17.61      3
#> 
#> $`0.0694444444444444`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00  33490
#> 2    0 0.00  31980
#> 3    0 0.00  30530
#> 4    0 0.00  29150
#> 5    0 0.00  27830
#> 6    0 0.00  26570
#> 7    0 0.00  25370
#> 8    0 0.00  24220
#> 9    0 0.00  23130
#> 10   0 0.00  22080
#> 11   0 0.00  21080
#> 12   0 0.00  20130
#> 13   0 0.00  19220
#> 14   0 0.00  18350
#> 15   0 0.00  17520
#> 16   0 0.00  16730
#> 17   0 0.00  15970
#> 18   0 0.00  15250
#> 19   0 0.00  14560
#> 20   0 0.00  13900
#> 21   0 0.00  13270
#> 22   0 0.00  12670
#> 23   0 0.00  12100
#> 24   0 0.00  11550
#> 25   0 0.00  11030
#> 26   0 0.00  10530
#> 27   0 0.00  10050
#> 28   0 0.00   9598
#> 29   0 0.00   9164
#> 30   0 0.00   8750
#> 31   0 0.00   8354
#> 32   0 0.00   7976
#> 33   0 0.00   7616
#> 34   0 0.00   7271
#> 35   0 0.00   6942
#> 36   0 0.00   6628
#> 37   0 0.00   6328
#> 38   0 0.00   6042
#> 39   0 0.00   5769
#> 40   0 0.00   5508
#> 41   0 0.00   5259
#> 42   0 0.00   5021
#> 43   0 0.00   4794
#> 44   0 0.00   4577
#> 45   0 0.00   4370
#> 46   0 0.00   4172
#> 47   0 0.00   3984
#> 48   0 0.00   3804
#> 49   0 0.00   3632
#> 50   0 0.00   3467
#> 51   0 0.00   3310
#> 52   0 0.00   3161
#> 53   0 0.00   3018
#> 54   0 0.00   2881
#> 55   0 0.00   2751
#> 56   0 0.00   2627
#> 57   0 0.00   2508
#> 58   0 0.00   2394
#> 59   0 0.00   2286
#> 60   0 0.00   2183
#> 61   0 0.00   2084
#> 62   0 0.00   1990
#> 63   0 0.00   1900
#> 64   0 0.00   1814
#> 65   0 0.00   1732
#> 66   0 0.00   1653
#> 67   0 0.00   1579
#> 68   0 0.00   1507
#> 69   0 0.00   1439
#> 70   0 0.00   1374
#> 71   0 0.00   1312
#> 72   0 0.00   1253
#> 73   0 0.00   1196
#> 74   0 0.00   1142
#> 75   0 0.00   1090
#> 76   0 0.00   1041
#> 77   0 0.00    994
#> 78   0 0.00    949
#> 79   0 0.00    906
#> 80   0 0.00    865
#> 81   0 0.00    826
#> 82   0 0.00    788
#> 83   0 0.00    753
#> 84   0 0.00    719
#> 85   0 0.00    686
#> 86   0 0.00    655
#> 87   0 0.00    626
#> 88   0 0.00    597
#> 89   0 0.00    570
#> 90   0 0.00    544
#> 91   0 0.00    520
#> 92   0 0.00    496
#> 93   0 0.00    474
#> 94   0 0.00    452
#> 95   0 0.00    432
#> 96   0 0.00    412
#> 97   0 0.00    394
#> 98   0 0.00    376
#> 99   0 0.00    359
#> 100  0 0.00    343
#> 101  0 0.00    327
#> 102  0 0.00    312
#> 103  0 0.00    298
#> 104  0 0.00    285
#> 105  0 0.00    272
#> 106  0 0.00    260
#> 107  0 0.00    248
#> 108  0 0.00    237
#> 109  0 0.00    226
#> 110  0 0.00    216
#> 111  0 0.00    206
#> 112  0 0.00    197
#> 113  0 0.00    188
#> 114  0 0.00    179
#> 115  0 0.00    171
#> 116  0 0.00    163
#> 117  0 0.00    156
#> 118  0 0.00    149
#> 119  0 0.00    142
#> 120  0 0.00    136
#> 121  0 0.00    130
#> 122  0 0.00    124
#> 123  0 0.00    118
#> 124  0 0.00    113
#> 125  0 0.00    108
#> 126  0 0.00    103
#> 127  0 0.00     98
#> 128  0 0.00     94
#> 129  0 0.00     90
#> 130  0 0.00     86
#> 131  0 0.00     82
#> 132  0 0.00     78
#> 133  0 0.00     74
#> 134  0 0.00     71
#> 135  0 0.00     68
#> 136  0 0.00     65
#> 137  0 0.00     62
#> 138  0 0.00     59
#> 139  0 0.00     56
#> 140  0 0.00     54
#> 141  0 0.00     51
#> 142  0 0.00     49
#> 143  0 0.00     47
#> 144  0 0.00     45
#> 145  0 0.00     43
#> 146  0 0.00     41
#> 147  0 0.00     39
#> 148  0 0.00     37
#> 149  0 0.00     35
#> 150  0 0.00     34
#> 151  0 0.00     32
#> 152  0 0.00     31
#> 153  0 0.00     29
#> 154  0 0.00     28
#> 155  0 0.00     27
#> 156  0 0.00     26
#> 157  0 0.00     24
#> 158  0 0.00     23
#> 159  0 0.00     22
#> 160  0 0.00     21
#> 161  0 0.00     20
#> 162  0 0.00     19
#> 163  0 0.00     19
#> 164  0 0.00     18
#> 165  0 0.00     17
#> 166  0 0.00     16
#> 167  0 0.00     15
#> 168  0 0.00     15
#> 169  0 0.00     14
#> 170  0 0.00     13
#> 171  0 0.00     13
#> 172  0 0.00     12
#> 173  0 0.00     12
#> 174  0 0.00     11
#> 175  0 0.00     11
#> 176  0 0.00     10
#> 177  1 0.09     10
#> 178  1 0.27      9
#> 179  1 0.43      9
#> 180  1 0.57      8
#> 181  1 0.70      8
#> 182  1 0.82      8
#> 183  1 0.92      7
#> 184  1 1.02      7
#> 185  1 1.10      7
#> 186  1 1.18      6
#> 187  1 1.25      6
#> 188  1 1.31      6
#> 189  1 1.37      6
#> 190  1 1.42      5
#> 191  1 1.47      5
#> 192  1 1.51      5
#> 193  1 1.55      5
#> 194  1 1.59      4
#> 195  1 1.62      4
#> 196  1 1.65      4
#> 197  1 1.67      4
#> 198  1 1.70      4
#> 199  1 1.72      4
#> 200  1 1.74      3
#> 
#> $`0.251157407407407`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00  33480
#> 2    0 0.00  31970
#> 3    0 0.00  30520
#> 4    0 0.00  29140
#> 5    0 0.00  27820
#> 6    0 0.00  26560
#> 7    0 0.00  25360
#> 8    0 0.00  24220
#> 9    0 0.00  23120
#> 10   0 0.00  22070
#> 11   0 0.00  21080
#> 12   0 0.00  20120
#> 13   0 0.00  19210
#> 14   0 0.00  18340
#> 15   0 0.00  17510
#> 16   0 0.00  16720
#> 17   0 0.00  15970
#> 18   0 0.00  15240
#> 19   0 0.00  14550
#> 20   0 0.00  13900
#> 21   0 0.00  13270
#> 22   0 0.00  12670
#> 23   0 0.00  12090
#> 24   0 0.00  11550
#> 25   0 0.00  11030
#> 26   0 0.00  10530
#> 27   0 0.00  10050
#> 28   0 0.00   9596
#> 29   0 0.00   9162
#> 30   0 0.00   8748
#> 31   0 0.00   8352
#> 32   0 0.00   7974
#> 33   0 0.00   7614
#> 34   0 0.00   7269
#> 35   0 0.00   6940
#> 36   0 0.00   6626
#> 37   0 0.00   6327
#> 38   0 0.00   6041
#> 39   0 0.00   5767
#> 40   0 0.00   5507
#> 41   0 0.00   5258
#> 42   0 0.00   5020
#> 43   0 0.00   4793
#> 44   0 0.00   4576
#> 45   0 0.00   4369
#> 46   0 0.00   4171
#> 47   0 0.00   3983
#> 48   0 0.00   3803
#> 49   0 0.00   3631
#> 50   0 0.00   3466
#> 51   0 0.00   3310
#> 52   0 0.00   3160
#> 53   0 0.00   3017
#> 54   0 0.00   2881
#> 55   0 0.00   2750
#> 56   0 0.00   2626
#> 57   0 0.00   2507
#> 58   0 0.00   2394
#> 59   0 0.00   2285
#> 60   0 0.00   2182
#> 61   0 0.00   2083
#> 62   0 0.00   1989
#> 63   0 0.00   1899
#> 64   0 0.00   1813
#> 65   0 0.00   1731
#> 66   0 0.00   1653
#> 67   0 0.00   1578
#> 68   0 0.00   1507
#> 69   0 0.00   1439
#> 70   0 0.00   1374
#> 71   0 0.00   1311
#> 72   0 0.00   1252
#> 73   0 0.00   1196
#> 74   0 0.00   1141
#> 75   0 0.00   1090
#> 76   0 0.00   1041
#> 77   0 0.00    994
#> 78   0 0.00    949
#> 79   0 0.00    906
#> 80   0 0.00    865
#> 81   0 0.00    826
#> 82   0 0.00    788
#> 83   0 0.00    753
#> 84   0 0.00    719
#> 85   0 0.00    686
#> 86   0 0.00    655
#> 87   0 0.00    625
#> 88   0 0.00    597
#> 89   0 0.00    570
#> 90   0 0.00    544
#> 91   0 0.00    520
#> 92   0 0.00    496
#> 93   0 0.00    474
#> 94   0 0.00    452
#> 95   0 0.00    432
#> 96   0 0.00    412
#> 97   0 0.00    394
#> 98   0 0.00    376
#> 99   0 0.00    359
#> 100  0 0.00    343
#> 101  0 0.00    327
#> 102  0 0.00    312
#> 103  0 0.00    298
#> 104  0 0.00    285
#> 105  0 0.00    272
#> 106  0 0.00    260
#> 107  0 0.00    248
#> 108  0 0.00    237
#> 109  0 0.00    226
#> 110  0 0.00    216
#> 111  0 0.00    206
#> 112  0 0.00    197
#> 113  0 0.00    188
#> 114  0 0.00    179
#> 115  0 0.00    171
#> 116  0 0.00    163
#> 117  0 0.00    156
#> 118  0 0.00    149
#> 119  0 0.00    142
#> 120  0 0.00    136
#> 121  0 0.00    130
#> 122  0 0.00    124
#> 123  0 0.00    118
#> 124  0 0.00    113
#> 125  0 0.00    108
#> 126  0 0.00    103
#> 127  0 0.00     98
#> 128  0 0.00     94
#> 129  0 0.00     90
#> 130  0 0.00     85
#> 131  0 0.00     82
#> 132  0 0.00     78
#> 133  0 0.00     74
#> 134  0 0.00     71
#> 135  0 0.00     68
#> 136  0 0.00     65
#> 137  0 0.00     62
#> 138  0 0.00     59
#> 139  0 0.00     56
#> 140  0 0.00     54
#> 141  0 0.00     51
#> 142  0 0.00     49
#> 143  0 0.00     47
#> 144  0 0.00     45
#> 145  0 0.00     43
#> 146  0 0.00     41
#> 147  0 0.00     39
#> 148  0 0.00     37
#> 149  0 0.00     35
#> 150  0 0.00     34
#> 151  0 0.00     32
#> 152  0 0.00     31
#> 153  0 0.00     29
#> 154  0 0.00     28
#> 155  0 0.00     27
#> 156  0 0.00     26
#> 157  0 0.00     24
#> 158  0 0.00     23
#> 159  0 0.00     22
#> 160  0 0.00     21
#> 161  1 0.13     20
#> 162  1 0.33     19
#> 163  1 0.53     19
#> 164  1 0.72     18
#> 165  1 0.91     17
#> 166  1 1.10     16
#> 167  1 1.28     15
#> 168  1 1.45     15
#> 169  1 1.62     14
#> 170  1 1.77     13
#> 171  1 1.92     13
#> 172  1 2.06     12
#> 173  1 2.20     12
#> 174  1 2.32     11
#> 175  1 2.44     11
#> 176  1 2.55     10
#> 177  1 2.65     10
#> 178  1 2.75      9
#> 179  1 2.83      9
#> 180  1 2.92      8
#> 181  1 2.99      8
#> 182  1 3.07      8
#> 183  1 3.13      7
#> 184  1 3.19      7
#> 185  1 3.25      7
#> 186  1 3.30      6
#> 187  1 3.35      6
#> 188  1 3.40      6
#> 189  1 3.44      6
#> 190  1 3.48      5
#> 191  1 3.52      5
#> 192  1 3.55      5
#> 193  1 3.58      5
#> 194  1 3.61      4
#> 195  1 3.63      4
#> 196  1 3.66      4
#> 197  1 3.68      4
#> 198  1 3.70      4
#> 199  1 3.72      4
#> 200  1 3.74      3
#> 
#> $`0.497685185185185`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00  42180
#> 2    0 0.00  40270
#> 3    0 0.00  38450
#> 4    0 0.00  36710
#> 5    0 0.00  35050
#> 6    0 0.00  33460
#> 7    0 0.00  31950
#> 8    0 0.00  30510
#> 9    0 0.00  29130
#> 10   0 0.00  27810
#> 11   0 0.00  26550
#> 12   0 0.00  25350
#> 13   0 0.00  24200
#> 14   0 0.00  23110
#> 15   0 0.00  22060
#> 16   0 0.00  21070
#> 17   0 0.00  20110
#> 18   0 0.00  19200
#> 19   0 0.00  18330
#> 20   0 0.00  17510
#> 21   0 0.00  16710
#> 22   0 0.00  15960
#> 23   0 0.00  15240
#> 24   0 0.00  14550
#> 25   0 0.00  13890
#> 26   0 0.00  13260
#> 27   0 0.00  12660
#> 28   0 0.00  12090
#> 29   0 0.00  11540
#> 30   0 0.00  11020
#> 31   0 0.00  10520
#> 32   0 0.00  10050
#> 33   0 0.00   9591
#> 34   0 0.00   9157
#> 35   0 0.00   8743
#> 36   0 0.00   8348
#> 37   0 0.00   7970
#> 38   0 0.00   7610
#> 39   0 0.00   7265
#> 40   0 0.00   6937
#> 41   0 0.00   6623
#> 42   0 0.00   6324
#> 43   0 0.00   6038
#> 44   0 0.00   5765
#> 45   0 0.00   5504
#> 46   0 0.00   5255
#> 47   0 0.00   5017
#> 48   0 0.00   4790
#> 49   0 0.00   4574
#> 50   0 0.00   4367
#> 51   0 0.00   4169
#> 52   0 0.00   3981
#> 53   0 0.00   3801
#> 54   0 0.00   3629
#> 55   0 0.00   3465
#> 56   0 0.00   3308
#> 57   0 0.00   3158
#> 58   0 0.00   3015
#> 59   0 0.00   2879
#> 60   0 0.00   2749
#> 61   0 0.00   2625
#> 62   0 0.00   2506
#> 63   0 0.00   2393
#> 64   0 0.00   2284
#> 65   0 0.00   2181
#> 66   0 0.00   2082
#> 67   0 0.00   1988
#> 68   0 0.00   1898
#> 69   0 0.00   1812
#> 70   0 0.00   1730
#> 71   0 0.00   1652
#> 72   0 0.00   1577
#> 73   0 0.00   1506
#> 74   0 0.00   1438
#> 75   0 0.00   1373
#> 76   0 0.00   1311
#> 77   0 0.00   1252
#> 78   0 0.00   1195
#> 79   0 0.00   1141
#> 80   0 0.00   1089
#> 81   0 0.00   1040
#> 82   0 0.00    993
#> 83   0 0.00    948
#> 84   0 0.00    905
#> 85   0 0.00    864
#> 86   0 0.00    825
#> 87   0 0.00    788
#> 88   0 0.00    752
#> 89   0 0.00    718
#> 90   0 0.00    686
#> 91   0 0.00    655
#> 92   0 0.00    625
#> 93   0 0.00    597
#> 94   0 0.00    570
#> 95   0 0.00    544
#> 96   0 0.00    519
#> 97   0 0.00    496
#> 98   0 0.00    474
#> 99   0 0.00    452
#> 100  0 0.00    432
#> 101  0 0.00    412
#> 102  0 0.00    394
#> 103  0 0.00    376
#> 104  0 0.00    359
#> 105  0 0.00    342
#> 106  0 0.00    327
#> 107  0 0.00    312
#> 108  0 0.00    298
#> 109  0 0.00    285
#> 110  0 0.00    272
#> 111  0 0.00    259
#> 112  0 0.00    248
#> 113  0 0.00    236
#> 114  0 0.00    226
#> 115  0 0.00    216
#> 116  0 0.00    206
#> 117  0 0.00    196
#> 118  0 0.00    188
#> 119  0 0.00    179
#> 120  0 0.00    171
#> 121  0 0.00    163
#> 122  0 0.00    156
#> 123  0 0.00    149
#> 124  0 0.00    142
#> 125  0 0.00    136
#> 126  0 0.00    130
#> 127  0 0.00    124
#> 128  0 0.00    118
#> 129  0 0.00    113
#> 130  0 0.00    108
#> 131  0 0.00    103
#> 132  0 0.00     98
#> 133  0 0.00     94
#> 134  0 0.00     89
#> 135  0 0.00     85
#> 136  0 0.00     82
#> 137  0 0.00     78
#> 138  0 0.00     74
#> 139  0 0.00     71
#> 140  0 0.00     68
#> 141  0 0.00     65
#> 142  0 0.00     62
#> 143  0 0.00     59
#> 144  0 0.00     56
#> 145  0 0.00     54
#> 146  0 0.00     51
#> 147  0 0.00     49
#> 148  0 0.00     47
#> 149  0 0.00     45
#> 150  0 0.00     43
#> 151  0 0.00     41
#> 152  0 0.00     39
#> 153  0 0.00     37
#> 154  0 0.00     35
#> 155  0 0.00     34
#> 156  0 0.00     32
#> 157  0 0.00     31
#> 158  0 0.00     29
#> 159  0 0.00     28
#> 160  0 0.00     27
#> 161  0 0.00     26
#> 162  0 0.00     24
#> 163  0 0.00     23
#> 164  0 0.00     22
#> 165  0 0.00     21
#> 166  0 0.00     20
#> 167  0 0.00     19
#> 168  0 0.00     19
#> 169  0 0.00     18
#> 170  0 0.00     17
#> 171  0 0.00     16
#> 172  0 0.00     15
#> 173  0 0.00     15
#> 174  0 0.00     14
#> 175  0 0.00     13
#> 176  0 0.00     13
#> 177  0 0.00     12
#> 178  0 0.00     12
#> 179  0 0.00     11
#> 180  0 0.00     11
#> 181  0 0.00     10
#> 182  0 0.00     10
#> 183  1 0.13      9
#> 184  1 0.27      9
#> 185  1 0.40      8
#> 186  1 0.52      8
#> 187  1 0.63      8
#> 188  1 0.73      7
#> 189  1 0.82      7
#> 190  1 0.90      7
#> 191  1 0.97      6
#> 192  1 1.04      6
#> 193  1 1.10      6
#> 194  1 1.16      6
#> 195  1 1.21      5
#> 196  1 1.26      5
#> 197  1 1.30      5
#> 198  1 1.34      5
#> 199  1 1.38      4
#> 200  1 1.41      4
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  64400
#> 2    0  0.00  61490
#> 3    0  0.00  58710
#> 4    0  0.00  56050
#> 5    0  0.00  53520
#> 6    0  0.00  51100
#> 7    0  0.00  48790
#> 8    0  0.00  46580
#> 9    0  0.00  44470
#> 10   0  0.00  42460
#> 11   0  0.00  40540
#> 12   0  0.00  38710
#> 13   0  0.00  36960
#> 14   0  0.00  35290
#> 15   0  0.00  33690
#> 16   0  0.00  32170
#> 17   0  0.00  30710
#> 18   0  0.00  29320
#> 19   0  0.00  28000
#> 20   0  0.00  26730
#> 21   0  0.00  25520
#> 22   0  0.00  24370
#> 23   0  0.00  23260
#> 24   0  0.00  22210
#> 25   0  0.00  21210
#> 26   0  0.00  20250
#> 27   0  0.00  19330
#> 28   0  0.00  18460
#> 29   0  0.00  17620
#> 30   0  0.00  16830
#> 31   0  0.00  16070
#> 32   0  0.00  15340
#> 33   0  0.00  14650
#> 34   0  0.00  13980
#> 35   0  0.00  13350
#> 36   0  0.00  12750
#> 37   0  0.00  12170
#> 38   0  0.00  11620
#> 39   0  0.00  11090
#> 40   0  0.00  10590
#> 41   0  0.00  10110
#> 42   0  0.00   9656
#> 43   0  0.00   9219
#> 44   0  0.00   8802
#> 45   0  0.00   8404
#> 46   0  0.00   8024
#> 47   0  0.00   7661
#> 48   0  0.00   7315
#> 49   0  0.00   6984
#> 50   0  0.00   6668
#> 51   0  0.00   6366
#> 52   0  0.00   6078
#> 53   0  0.00   5803
#> 54   0  0.00   5541
#> 55   0  0.00   5290
#> 56   0  0.00   5051
#> 57   0  0.00   4823
#> 58   0  0.00   4604
#> 59   0  0.00   4396
#> 60   0  0.00   4197
#> 61   0  0.00   4008
#> 62   0  0.00   3826
#> 63   0  0.00   3653
#> 64   0  0.00   3488
#> 65   0  0.00   3330
#> 66   0  0.00   3180
#> 67   0  0.00   3036
#> 68   0  0.00   2899
#> 69   0  0.00   2767
#> 70   0  0.00   2642
#> 71   0  0.00   2523
#> 72   0  0.00   2409
#> 73   0  0.00   2300
#> 74   0  0.00   2196
#> 75   0  0.00   2096
#> 76   0  0.00   2002
#> 77   0  0.00   1911
#> 78   0  0.00   1825
#> 79   0  0.00   1742
#> 80   0  0.00   1663
#> 81   0  0.00   1588
#> 82   0  0.00   1516
#> 83   0  0.00   1448
#> 84   0  0.00   1382
#> 85   0  0.00   1320
#> 86   0  0.00   1260
#> 87   0  0.00   1203
#> 88   0  0.00   1149
#> 89   0  0.00   1097
#> 90   0  0.00   1047
#> 91   0  0.00   1000
#> 92   0  0.00    954
#> 93   0  0.00    911
#> 94   0  0.00    870
#> 95   0  0.00    831
#> 96   0  0.00    793
#> 97   0  0.00    757
#> 98   0  0.00    723
#> 99   0  0.00    690
#> 100  0  0.00    659
#> 101  0  0.00    629
#> 102  0  0.00    601
#> 103  0  0.00    574
#> 104  0  0.00    548
#> 105  1  0.55    523
#> 106  1  1.35    499
#> 107  1  2.09    477
#> 108  1  2.78    455
#> 109  1  3.42    435
#> 110  1  4.02    415
#> 111  1  4.59    396
#> 112  1  5.12    378
#> 113  1  5.62    361
#> 114  1  6.10    345
#> 115  1  6.54    329
#> 116  1  6.96    314
#> 117  1  7.36    300
#> 118  1  7.73    286
#> 119  1  8.09    274
#> 120  1  8.43    261
#> 121  1  8.76    249
#> 122  1  9.07    238
#> 123  1  9.36    227
#> 124  1  9.64    217
#> 125  1  9.91    207
#> 126  1 10.17    198
#> 127  1 10.41    189
#> 128  1 10.65    180
#> 129  1 10.87    172
#> 130  1 11.08    164
#> 131  1 11.29    157
#> 132  1 11.49    150
#> 133  1 11.68    143
#> 134  1 11.86    137
#> 135  1 12.03    130
#> 136  1 12.20    124
#> 137  1 12.36    119
#> 138  1 12.51    114
#> 139  1 12.66    108
#> 140  1 12.80    104
#> 141  1 12.94     99
#> 142  1 13.07     94
#> 143  1 13.20     90
#> 144  1 13.32     86
#> 145  1 13.44     82
#> 146  1 13.55     78
#> 147  1 13.66     75
#> 148  1 13.76     71
#> 149  1 13.86     68
#> 150  1 13.96     65
#> 151  1 14.05     62
#> 152  1 14.14     59
#> 153  1 14.22     57
#> 154  1 14.31     54
#> 155  1 14.38     52
#> 156  1 14.46     49
#> 157  1 14.53     47
#> 158  1 14.60     45
#> 159  1 14.66     43
#> 160  1 14.73     41
#> 161  1 14.79     39
#> 162  1 14.85     37
#> 163  1 14.90     36
#> 164  1 14.95     34
#> 165  1 15.00     33
#> 166  1 15.05     31
#> 167  1 15.09     30
#> 168  1 15.14     28
#> 169  1 15.18     27
#> 170  1 15.22     26
#> 171  1 15.25     25
#> 172  1 15.29     24
#> 173  1 15.32     22
#> 174  1 15.35     21
#> 175  1 15.38     20
#> 176  1 15.41     20
#> 177  1 15.43     19
#> 178  1 15.46     18
#> 179  1 15.48     17
#> 180  1 15.50     16
#> 181  2 15.54     16
#> 182  2 15.63     15
#> 183  2 15.71     14
#> 184  2 15.79     14
#> 185  2 15.87     13
#> 186  2 15.93     12
#> 187  2 16.00     12
#> 188  2 16.06     11
#> 189  2 16.12     11
#> 190  2 16.17     10
#> 191  2 16.22     10
#> 192  2 16.26      9
#> 193  2 16.31      9
#> 194  2 16.35      9
#> 195  2 16.38      8
#> 196  2 16.42      8
#> 197  2 16.45      7
#> 198  2 16.48      7
#> 199  2 16.51      7
#> 200  2 16.53      6
#> 
#> $`0.0578703703703702`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  90880
#> 2    0  0.00  86770
#> 3    0  0.00  82840
#> 4    0  0.00  79100
#> 5    0  0.00  75520
#> 6    0  0.00  72100
#> 7    0  0.00  68840
#> 8    0  0.00  65730
#> 9    0  0.00  62760
#> 10   0  0.00  59920
#> 11   0  0.00  57210
#> 12   0  0.00  54620
#> 13   0  0.00  52150
#> 14   0  0.00  49790
#> 15   0  0.00  47540
#> 16   0  0.00  45390
#> 17   0  0.00  43340
#> 18   0  0.00  41380
#> 19   0  0.00  39500
#> 20   0  0.00  37720
#> 21   0  0.00  36010
#> 22   0  0.00  34380
#> 23   0  0.00  32830
#> 24   0  0.00  31340
#> 25   0  0.00  29930
#> 26   0  0.00  28570
#> 27   0  0.00  27280
#> 28   0  0.00  26050
#> 29   0  0.00  24870
#> 30   0  0.00  23740
#> 31   0  0.00  22670
#> 32   0  0.00  21640
#> 33   0  0.00  20670
#> 34   0  0.00  19730
#> 35   0  0.00  18840
#> 36   0  0.00  17990
#> 37   0  0.00  17170
#> 38   0  0.00  16400
#> 39   0  0.00  15650
#> 40   0  0.00  14950
#> 41   0  0.00  14270
#> 42   0  0.00  13620
#> 43   0  0.00  13010
#> 44   0  0.00  12420
#> 45   0  0.00  11860
#> 46   0  0.00  11320
#> 47   0  0.00  10810
#> 48   0  0.00  10320
#> 49   0  0.00   9854
#> 50   0  0.00   9409
#> 51   0  0.00   8983
#> 52   0  0.00   8577
#> 53   0  0.00   8189
#> 54   0  0.00   7819
#> 55   0  0.00   7465
#> 56   0  0.00   7127
#> 57   0  0.00   6805
#> 58   0  0.00   6497
#> 59   0  0.00   6203
#> 60   0  0.00   5923
#> 61   0  0.00   5655
#> 62   0  0.00   5399
#> 63   0  0.00   5155
#> 64   0  0.00   4922
#> 65   0  0.00   4699
#> 66   0  0.00   4487
#> 67   0  0.00   4284
#> 68   0  0.00   4090
#> 69   0  0.00   3905
#> 70   0  0.00   3728
#> 71   0  0.00   3560
#> 72   0  0.00   3399
#> 73   0  0.00   3245
#> 74   0  0.00   3098
#> 75   0  0.00   2958
#> 76   0  0.00   2824
#> 77   0  0.00   2697
#> 78   0  0.00   2575
#> 79   0  0.00   2458
#> 80   0  0.00   2347
#> 81   0  0.00   2241
#> 82   0  0.00   2140
#> 83   0  0.00   2043
#> 84   0  0.00   1950
#> 85   0  0.00   1862
#> 86   0  0.00   1778
#> 87   0  0.00   1698
#> 88   0  0.00   1621
#> 89   0  0.00   1547
#> 90   0  0.00   1477
#> 91   0  0.00   1411
#> 92   0  0.00   1347
#> 93   0  0.00   1286
#> 94   0  0.00   1228
#> 95   0  0.00   1172
#> 96   0  0.00   1119
#> 97   0  0.00   1069
#> 98   0  0.00   1020
#> 99   0  0.00    974
#> 100  0  0.00    930
#> 101  0  0.00    888
#> 102  0  0.00    848
#> 103  0  0.00    810
#> 104  0  0.00    773
#> 105  0  0.00    738
#> 106  0  0.00    704
#> 107  0  0.00    673
#> 108  0  0.00    642
#> 109  0  0.00    613
#> 110  0  0.00    586
#> 111  0  0.00    559
#> 112  0  0.00    534
#> 113  0  0.00    510
#> 114  0  0.00    486
#> 115  0  0.00    464
#> 116  0  0.00    444
#> 117  0  0.00    423
#> 118  0  0.00    404
#> 119  0  0.00    386
#> 120  0  0.00    368
#> 121  0  0.00    352
#> 122  0  0.00    336
#> 123  0  0.00    321
#> 124  0  0.00    306
#> 125  0  0.00    292
#> 126  0  0.00    279
#> 127  0  0.00    267
#> 128  0  0.00    254
#> 129  0  0.00    243
#> 130  0  0.00    232
#> 131  0  0.00    222
#> 132  0  0.00    212
#> 133  0  0.00    202
#> 134  0  0.00    193
#> 135  0  0.00    184
#> 136  0  0.00    176
#> 137  0  0.00    168
#> 138  0  0.00    160
#> 139  0  0.00    153
#> 140  1  0.10    146
#> 141  1  0.77    139
#> 142  1  1.39    133
#> 143  1  1.97    127
#> 144  1  2.52    121
#> 145  1  3.03    116
#> 146  1  3.51    111
#> 147  1  3.97    106
#> 148  1  4.39    101
#> 149  1  4.80     96
#> 150  1  5.18     92
#> 151  1  5.54     88
#> 152  1  5.87     84
#> 153  1  6.20     80
#> 154  1  6.50     76
#> 155  1  6.79     73
#> 156  1  7.06     70
#> 157  1  7.32     66
#> 158  1  7.56     63
#> 159  1  7.79     61
#> 160  1  8.01     58
#> 161  1  8.21     55
#> 162  1  8.41     53
#> 163  1  8.59     50
#> 164  1  8.77     48
#> 165  1  8.93     46
#> 166  1  9.09     44
#> 167  1  9.24     42
#> 168  1  9.38     40
#> 169  1  9.51     38
#> 170  1  9.64     36
#> 171  1  9.75     35
#> 172  1  9.86     33
#> 173  1  9.97     32
#> 174  1 10.06     30
#> 175  1 10.16     29
#> 176  1 10.24     28
#> 177  1 10.33     26
#> 178  1 10.40     25
#> 179  1 10.47     24
#> 180  1 10.54     23
#> 181  1 10.61     22
#> 182  1 10.66     21
#> 183  1 10.72     20
#> 184  1 10.77     19
#> 185  1 10.82     18
#> 186  1 10.86     17
#> 187  1 10.91     17
#> 188  1 10.95     16
#> 189  1 10.98     15
#> 190  1 11.02     14
#> 191  1 11.05     14
#> 192  1 11.08     13
#> 193  1 11.10     13
#> 194  1 11.13     12
#> 195  1 11.15     11
#> 196  1 11.17     11
#> 197  1 11.19     10
#> 198  1 11.21     10
#> 199  1 11.23     10
#> 200  1 11.25      9
#> 
#> $`0.0578703703703702`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  64520
#> 2    0  0.00  61600
#> 3    0  0.00  58810
#> 4    0  0.00  56150
#> 5    0  0.00  53610
#> 6    0  0.00  51190
#> 7    0  0.00  48870
#> 8    0  0.00  46660
#> 9    0  0.00  44550
#> 10   0  0.00  42540
#> 11   0  0.00  40610
#> 12   0  0.00  38780
#> 13   0  0.00  37020
#> 14   0  0.00  35350
#> 15   0  0.00  33750
#> 16   0  0.00  32220
#> 17   0  0.00  30770
#> 18   0  0.00  29370
#> 19   0  0.00  28050
#> 20   0  0.00  26780
#> 21   0  0.00  25570
#> 22   0  0.00  24410
#> 23   0  0.00  23310
#> 24   0  0.00  22250
#> 25   0  0.00  21250
#> 26   0  0.00  20280
#> 27   0  0.00  19370
#> 28   0  0.00  18490
#> 29   0  0.00  17650
#> 30   0  0.00  16860
#> 31   0  0.00  16090
#> 32   0  0.00  15370
#> 33   0  0.00  14670
#> 34   0  0.00  14010
#> 35   0  0.00  13370
#> 36   0  0.00  12770
#> 37   0  0.00  12190
#> 38   0  0.00  11640
#> 39   0  0.00  11110
#> 40   0  0.00  10610
#> 41   0  0.00  10130
#> 42   0  0.00   9673
#> 43   0  0.00   9235
#> 44   0  0.00   8818
#> 45   0  0.00   8419
#> 46   0  0.00   8038
#> 47   0  0.00   7674
#> 48   0  0.00   7327
#> 49   0  0.00   6996
#> 50   0  0.00   6680
#> 51   0  0.00   6377
#> 52   0  0.00   6089
#> 53   0  0.00   5814
#> 54   0  0.00   5551
#> 55   0  0.00   5300
#> 56   0  0.00   5060
#> 57   0  0.00   4831
#> 58   0  0.00   4613
#> 59   0  0.00   4404
#> 60   0  0.00   4205
#> 61   0  0.00   4015
#> 62   0  0.00   3833
#> 63   0  0.00   3660
#> 64   0  0.00   3494
#> 65   0  0.00   3336
#> 66   0  0.00   3185
#> 67   0  0.00   3041
#> 68   0  0.00   2904
#> 69   0  0.00   2772
#> 70   0  0.00   2647
#> 71   0  0.00   2527
#> 72   0  0.00   2413
#> 73   0  0.00   2304
#> 74   0  0.00   2200
#> 75   0  0.00   2100
#> 76   0  0.00   2005
#> 77   0  0.00   1914
#> 78   0  0.00   1828
#> 79   0  0.00   1745
#> 80   0  0.00   1666
#> 81   0  0.00   1591
#> 82   0  0.00   1519
#> 83   0  0.00   1450
#> 84   0  0.00   1385
#> 85   0  0.00   1322
#> 86   0  0.00   1262
#> 87   0  0.00   1205
#> 88   0  0.00   1151
#> 89   0  0.00   1099
#> 90   0  0.00   1049
#> 91   0  0.00   1001
#> 92   0  0.00    956
#> 93   0  0.00    913
#> 94   0  0.00    872
#> 95   0  0.00    832
#> 96   0  0.00    795
#> 97   0  0.00    759
#> 98   0  0.00    724
#> 99   0  0.00    692
#> 100  0  0.00    660
#> 101  0  0.00    630
#> 102  0  0.00    602
#> 103  0  0.00    575
#> 104  0  0.00    549
#> 105  0  0.00    524
#> 106  0  0.00    500
#> 107  0  0.00    478
#> 108  0  0.00    456
#> 109  0  0.00    435
#> 110  0  0.00    416
#> 111  0  0.00    397
#> 112  0  0.00    379
#> 113  0  0.00    362
#> 114  0  0.00    345
#> 115  0  0.00    330
#> 116  0  0.00    315
#> 117  0  0.00    301
#> 118  0  0.00    287
#> 119  0  0.00    274
#> 120  0  0.00    262
#> 121  0  0.00    250
#> 122  0  0.00    238
#> 123  0  0.00    228
#> 124  0  0.00    217
#> 125  0  0.00    208
#> 126  0  0.00    198
#> 127  0  0.00    189
#> 128  0  0.00    181
#> 129  0  0.00    172
#> 130  0  0.00    165
#> 131  0  0.00    157
#> 132  0  0.00    150
#> 133  0  0.00    143
#> 134  0  0.00    137
#> 135  0  0.00    131
#> 136  0  0.00    125
#> 137  0  0.00    119
#> 138  0  0.00    114
#> 139  0  0.00    109
#> 140  1  0.10    104
#> 141  1  0.77     99
#> 142  1  1.39     95
#> 143  1  1.97     90
#> 144  1  2.52     86
#> 145  1  3.03     82
#> 146  1  3.51     79
#> 147  1  3.97     75
#> 148  1  4.39     72
#> 149  1  4.80     68
#> 150  1  5.18     65
#> 151  1  5.54     62
#> 152  1  5.87     60
#> 153  1  6.20     57
#> 154  1  6.50     54
#> 155  1  6.79     52
#> 156  1  7.06     49
#> 157  1  7.32     47
#> 158  1  7.56     45
#> 159  1  7.79     43
#> 160  1  8.01     41
#> 161  1  8.21     39
#> 162  1  8.41     37
#> 163  1  8.59     36
#> 164  1  8.77     34
#> 165  1  8.93     33
#> 166  1  9.09     31
#> 167  1  9.24     30
#> 168  1  9.38     28
#> 169  1  9.51     27
#> 170  1  9.64     26
#> 171  1  9.75     25
#> 172  1  9.86     24
#> 173  1  9.97     23
#> 174  1 10.06     21
#> 175  1 10.16     21
#> 176  1 10.24     20
#> 177  1 10.33     19
#> 178  1 10.40     18
#> 179  1 10.47     17
#> 180  1 10.54     16
#> 181  1 10.61     16
#> 182  1 10.66     15
#> 183  1 10.72     14
#> 184  1 10.77     14
#> 185  1 10.82     13
#> 186  1 10.86     12
#> 187  1 10.91     12
#> 188  1 10.95     11
#> 189  1 10.98     11
#> 190  1 11.02     10
#> 191  1 11.05     10
#> 192  1 11.08      9
#> 193  1 11.10      9
#> 194  1 11.13      9
#> 195  1 11.15      8
#> 196  1 11.17      8
#> 197  1 11.19      7
#> 198  1 11.21      7
#> 199  1 11.23      7
#> 200  1 11.25      6
#> 
#> $`0.0208333333333333`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  64410
#> 2    0  0.00  61500
#> 3    0  0.00  58720
#> 4    0  0.00  56060
#> 5    0  0.00  53520
#> 6    0  0.00  51100
#> 7    0  0.00  48790
#> 8    0  0.00  46590
#> 9    0  0.00  44480
#> 10   0  0.00  42470
#> 11   0  0.00  40550
#> 12   0  0.00  38710
#> 13   0  0.00  36960
#> 14   0  0.00  35290
#> 15   0  0.00  33690
#> 16   0  0.00  32170
#> 17   0  0.00  30710
#> 18   0  0.00  29330
#> 19   0  0.00  28000
#> 20   0  0.00  26730
#> 21   0  0.00  25520
#> 22   0  0.00  24370
#> 23   0  0.00  23270
#> 24   0  0.00  22210
#> 25   0  0.00  21210
#> 26   0  0.00  20250
#> 27   0  0.00  19330
#> 28   0  0.00  18460
#> 29   0  0.00  17630
#> 30   0  0.00  16830
#> 31   0  0.00  16070
#> 32   0  0.00  15340
#> 33   0  0.00  14650
#> 34   0  0.00  13980
#> 35   0  0.00  13350
#> 36   0  0.00  12750
#> 37   0  0.00  12170
#> 38   0  0.00  11620
#> 39   0  0.00  11100
#> 40   0  0.00  10590
#> 41   0  0.00  10110
#> 42   0  0.00   9657
#> 43   0  0.00   9220
#> 44   0  0.00   8803
#> 45   0  0.00   8405
#> 46   0  0.00   8025
#> 47   0  0.00   7662
#> 48   0  0.00   7315
#> 49   0  0.00   6984
#> 50   0  0.00   6669
#> 51   0  0.00   6367
#> 52   0  0.00   6079
#> 53   0  0.00   5804
#> 54   0  0.00   5542
#> 55   0  0.00   5291
#> 56   0  0.00   5052
#> 57   0  0.00   4823
#> 58   0  0.00   4605
#> 59   0  0.00   4397
#> 60   0  0.00   4198
#> 61   0  0.00   4008
#> 62   0  0.00   3827
#> 63   0  0.00   3654
#> 64   0  0.00   3488
#> 65   0  0.00   3331
#> 66   0  0.00   3180
#> 67   0  0.00   3036
#> 68   0  0.00   2899
#> 69   0  0.00   2768
#> 70   0  0.00   2643
#> 71   0  0.00   2523
#> 72   0  0.00   2409
#> 73   0  0.00   2300
#> 74   0  0.00   2196
#> 75   0  0.00   2097
#> 76   0  0.00   2002
#> 77   0  0.00   1911
#> 78   0  0.00   1825
#> 79   0  0.00   1742
#> 80   0  0.00   1663
#> 81   0  0.00   1588
#> 82   0  0.00   1516
#> 83   0  0.00   1448
#> 84   0  0.00   1382
#> 85   0  0.00   1320
#> 86   0  0.00   1260
#> 87   0  0.00   1203
#> 88   0  0.00   1149
#> 89   0  0.00   1097
#> 90   0  0.00   1047
#> 91   0  0.00   1000
#> 92   0  0.00    955
#> 93   0  0.00    911
#> 94   0  0.00    870
#> 95   0  0.00    831
#> 96   0  0.00    793
#> 97   0  0.00    757
#> 98   0  0.00    723
#> 99   0  0.00    690
#> 100  0  0.00    659
#> 101  0  0.00    629
#> 102  0  0.00    601
#> 103  0  0.00    574
#> 104  0  0.00    548
#> 105  0  0.00    523
#> 106  0  0.00    499
#> 107  0  0.00    477
#> 108  0  0.00    455
#> 109  0  0.00    435
#> 110  0  0.00    415
#> 111  0  0.00    396
#> 112  0  0.00    378
#> 113  0  0.00    361
#> 114  0  0.00    345
#> 115  0  0.00    329
#> 116  0  0.00    314
#> 117  0  0.00    300
#> 118  0  0.00    286
#> 119  0  0.00    274
#> 120  0  0.00    261
#> 121  0  0.00    249
#> 122  0  0.00    238
#> 123  0  0.00    227
#> 124  0  0.00    217
#> 125  0  0.00    207
#> 126  0  0.00    198
#> 127  0  0.00    189
#> 128  0  0.00    180
#> 129  0  0.00    172
#> 130  0  0.00    164
#> 131  0  0.00    157
#> 132  0  0.00    150
#> 133  0  0.00    143
#> 134  0  0.00    137
#> 135  0  0.00    130
#> 136  0  0.00    125
#> 137  0  0.00    119
#> 138  0  0.00    114
#> 139  0  0.00    108
#> 140  1  0.10    104
#> 141  1  0.77     99
#> 142  1  1.39     94
#> 143  1  1.97     90
#> 144  1  2.52     86
#> 145  1  3.03     82
#> 146  1  3.51     78
#> 147  1  3.97     75
#> 148  1  4.39     71
#> 149  1  4.80     68
#> 150  1  5.18     65
#> 151  1  5.54     62
#> 152  1  5.87     59
#> 153  1  6.20     57
#> 154  1  6.50     54
#> 155  1  6.79     52
#> 156  1  7.06     49
#> 157  1  7.32     47
#> 158  1  7.56     45
#> 159  1  7.79     43
#> 160  1  8.01     41
#> 161  1  8.21     39
#> 162  1  8.41     37
#> 163  1  8.59     36
#> 164  1  8.77     34
#> 165  1  8.93     33
#> 166  1  9.09     31
#> 167  1  9.24     30
#> 168  1  9.38     28
#> 169  1  9.51     27
#> 170  1  9.64     26
#> 171  1  9.75     25
#> 172  1  9.86     24
#> 173  1  9.97     22
#> 174  1 10.06     21
#> 175  1 10.16     20
#> 176  1 10.24     20
#> 177  1 10.33     19
#> 178  1 10.40     18
#> 179  1 10.47     17
#> 180  1 10.54     16
#> 181  1 10.61     16
#> 182  1 10.66     15
#> 183  1 10.72     14
#> 184  1 10.77     14
#> 185  1 10.82     13
#> 186  1 10.86     12
#> 187  1 10.91     12
#> 188  1 10.95     11
#> 189  1 10.98     11
#> 190  1 11.02     10
#> 191  1 11.05     10
#> 192  1 11.08      9
#> 193  1 11.10      9
#> 194  1 11.13      9
#> 195  1 11.15      8
#> 196  1 11.17      8
#> 197  1 11.19      7
#> 198  1 11.21      7
#> 199  1 11.23      7
#> 200  1 11.25      6
#> 
#> $`0.141203703703704`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  64400
#> 2    0  0.00  61490
#> 3    0  0.00  58710
#> 4    0  0.00  56050
#> 5    0  0.00  53520
#> 6    0  0.00  51100
#> 7    0  0.00  48790
#> 8    0  0.00  46580
#> 9    0  0.00  44470
#> 10   0  0.00  42460
#> 11   0  0.00  40540
#> 12   0  0.00  38710
#> 13   0  0.00  36960
#> 14   0  0.00  35280
#> 15   0  0.00  33690
#> 16   0  0.00  32170
#> 17   0  0.00  30710
#> 18   0  0.00  29320
#> 19   0  0.00  28000
#> 20   0  0.00  26730
#> 21   0  0.00  25520
#> 22   0  0.00  24370
#> 23   0  0.00  23260
#> 24   0  0.00  22210
#> 25   0  0.00  21210
#> 26   0  0.00  20250
#> 27   0  0.00  19330
#> 28   0  0.00  18460
#> 29   0  0.00  17620
#> 30   0  0.00  16830
#> 31   0  0.00  16060
#> 32   0  0.00  15340
#> 33   0  0.00  14640
#> 34   0  0.00  13980
#> 35   0  0.00  13350
#> 36   0  0.00  12750
#> 37   0  0.00  12170
#> 38   0  0.00  11620
#> 39   0  0.00  11090
#> 40   0  0.00  10590
#> 41   0  0.00  10110
#> 42   0  0.00   9655
#> 43   0  0.00   9219
#> 44   0  0.00   8802
#> 45   0  0.00   8404
#> 46   0  0.00   8024
#> 47   0  0.00   7661
#> 48   0  0.00   7314
#> 49   0  0.00   6983
#> 50   0  0.00   6668
#> 51   0  0.00   6366
#> 52   0  0.00   6078
#> 53   0  0.00   5803
#> 54   0  0.00   5541
#> 55   0  0.00   5290
#> 56   0  0.00   5051
#> 57   0  0.00   4822
#> 58   0  0.00   4604
#> 59   0  0.00   4396
#> 60   0  0.00   4197
#> 61   0  0.00   4007
#> 62   0  0.00   3826
#> 63   0  0.00   3653
#> 64   0  0.00   3488
#> 65   0  0.00   3330
#> 66   0  0.00   3180
#> 67   0  0.00   3036
#> 68   0  0.00   2898
#> 69   0  0.00   2767
#> 70   0  0.00   2642
#> 71   0  0.00   2523
#> 72   0  0.00   2409
#> 73   0  0.00   2300
#> 74   0  0.00   2196
#> 75   0  0.00   2096
#> 76   0  0.00   2002
#> 77   0  0.00   1911
#> 78   0  0.00   1825
#> 79   0  0.00   1742
#> 80   0  0.00   1663
#> 81   0  0.00   1588
#> 82   0  0.00   1516
#> 83   0  0.00   1448
#> 84   0  0.00   1382
#> 85   0  0.00   1320
#> 86   0  0.00   1260
#> 87   0  0.00   1203
#> 88   0  0.00   1149
#> 89   0  0.00   1097
#> 90   0  0.00   1047
#> 91   0  0.00   1000
#> 92   0  0.00    954
#> 93   0  0.00    911
#> 94   0  0.00    870
#> 95   0  0.00    831
#> 96   0  0.00    793
#> 97   0  0.00    757
#> 98   0  0.00    723
#> 99   0  0.00    690
#> 100  0  0.00    659
#> 101  0  0.00    629
#> 102  0  0.00    601
#> 103  0  0.00    574
#> 104  0  0.00    548
#> 105  0  0.00    523
#> 106  0  0.00    499
#> 107  0  0.00    477
#> 108  0  0.00    455
#> 109  0  0.00    435
#> 110  0  0.00    415
#> 111  0  0.00    396
#> 112  0  0.00    378
#> 113  0  0.00    361
#> 114  0  0.00    345
#> 115  0  0.00    329
#> 116  0  0.00    314
#> 117  0  0.00    300
#> 118  0  0.00    286
#> 119  0  0.00    274
#> 120  0  0.00    261
#> 121  0  0.00    249
#> 122  0  0.00    238
#> 123  0  0.00    227
#> 124  0  0.00    217
#> 125  0  0.00    207
#> 126  0  0.00    198
#> 127  0  0.00    189
#> 128  0  0.00    180
#> 129  0  0.00    172
#> 130  0  0.00    164
#> 131  0  0.00    157
#> 132  0  0.00    150
#> 133  0  0.00    143
#> 134  0  0.00    137
#> 135  0  0.00    130
#> 136  0  0.00    124
#> 137  0  0.00    119
#> 138  0  0.00    114
#> 139  0  0.00    108
#> 140  1  0.10    104
#> 141  1  0.77     99
#> 142  1  1.39     94
#> 143  1  1.97     90
#> 144  1  2.52     86
#> 145  1  3.03     82
#> 146  1  3.51     78
#> 147  1  3.97     75
#> 148  1  4.39     71
#> 149  1  4.80     68
#> 150  1  5.18     65
#> 151  1  5.54     62
#> 152  1  5.87     59
#> 153  1  6.20     57
#> 154  1  6.50     54
#> 155  1  6.79     52
#> 156  1  7.06     49
#> 157  1  7.32     47
#> 158  1  7.56     45
#> 159  1  7.79     43
#> 160  1  8.01     41
#> 161  1  8.21     39
#> 162  1  8.41     37
#> 163  1  8.59     36
#> 164  1  8.77     34
#> 165  1  8.93     33
#> 166  1  9.09     31
#> 167  1  9.24     30
#> 168  1  9.38     28
#> 169  1  9.51     27
#> 170  1  9.64     26
#> 171  1  9.75     25
#> 172  1  9.86     24
#> 173  1  9.97     22
#> 174  1 10.06     21
#> 175  1 10.16     20
#> 176  1 10.24     20
#> 177  1 10.33     19
#> 178  1 10.40     18
#> 179  1 10.47     17
#> 180  2 10.55     16
#> 181  2 10.64     16
#> 182  2 10.73     15
#> 183  2 10.81     14
#> 184  2 10.88     14
#> 185  2 10.95     13
#> 186  2 11.01     12
#> 187  2 11.08     12
#> 188  2 11.13     11
#> 189  2 11.18     11
#> 190  2 11.23     10
#> 191  2 11.28     10
#> 192  2 11.32      9
#> 193  2 11.36      9
#> 194  2 11.40      9
#> 195  2 11.43      8
#> 196  2 11.46      8
#> 197  2 11.49      7
#> 198  2 11.52      7
#> 199  2 11.54      7
#> 200  2 11.57      6
#> 
#> $`0.141203703703704`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00  73100
#> 2    0  0.00  69790
#> 3    0  0.00  66640
#> 4    0  0.00  63620
#> 5    0  0.00  60740
#> 6    0  0.00  58000
#> 7    0  0.00  55370
#> 8    0  0.00  52870
#> 9    0  0.00  50480
#> 10   0  0.00  48190
#> 11   0  0.00  46020
#> 12   0  0.00  43930
#> 13   0  0.00  41950
#> 14   0  0.00  40050
#> 15   0  0.00  38240
#> 16   0  0.00  36510
#> 17   0  0.00  34860
#> 18   0  0.00  33280
#> 19   0  0.00  31780
#> 20   0  0.00  30340
#> 21   0  0.00  28970
#> 22   0  0.00  27660
#> 23   0  0.00  26410
#> 24   0  0.00  25210
#> 25   0  0.00  24070
#> 26   0  0.00  22980
#> 27   0  0.00  21940
#> 28   0  0.00  20950
#> 29   0  0.00  20000
#> 30   0  0.00  19100
#> 31   0  0.00  18230
#> 32   0  0.00  17410
#> 33   0  0.00  16620
#> 34   0  0.00  15870
#> 35   0  0.00  15150
#> 36   0  0.00  14470
#> 37   0  0.00  13810
#> 38   0  0.00  13190
#> 39   0  0.00  12590
#> 40   0  0.00  12020
#> 41   0  0.00  11480
#> 42   0  0.00  10960
#> 43   0  0.00  10460
#> 44   0  0.00   9990
#> 45   0  0.00   9539
#> 46   0  0.00   9107
#> 47   0  0.00   8695
#> 48   0  0.00   8302
#> 49   0  0.00   7926
#> 50   0  0.00   7568
#> 51   0  0.00   7226
#> 52   0  0.00   6899
#> 53   0  0.00   6587
#> 54   0  0.00   6289
#> 55   0  0.00   6005
#> 56   0  0.00   5733
#> 57   0  0.00   5474
#> 58   0  0.00   5226
#> 59   0  0.00   4990
#> 60   0  0.00   4764
#> 61   0  0.00   4549
#> 62   0  0.00   4343
#> 63   0  0.00   4146
#> 64   0  0.00   3959
#> 65   0  0.00   3780
#> 66   0  0.00   3609
#> 67   0  0.00   3446
#> 68   0  0.00   3290
#> 69   0  0.00   3141
#> 70   0  0.00   2999
#> 71   0  0.00   2863
#> 72   0  0.00   2734
#> 73   0  0.00   2610
#> 74   0  0.00   2492
#> 75   0  0.00   2379
#> 76   0  0.00   2272
#> 77   0  0.00   2169
#> 78   0  0.00   2071
#> 79   0  0.00   1977
#> 80   0  0.00   1888
#> 81   0  0.00   1802
#> 82   0  0.00   1721
#> 83   0  0.00   1643
#> 84   0  0.00   1569
#> 85   0  0.00   1498
#> 86   0  0.00   1430
#> 87   0  0.00   1365
#> 88   0  0.00   1304
#> 89   0  0.00   1245
#> 90   0  0.00   1188
#> 91   0  0.00   1135
#> 92   0  0.00   1083
#> 93   0  0.00   1034
#> 94   0  0.00    988
#> 95   0  0.00    943
#> 96   0  0.00    900
#> 97   0  0.00    860
#> 98   0  0.00    821
#> 99   0  0.00    784
#> 100  0  0.00    748
#> 101  0  0.00    714
#> 102  0  0.00    682
#> 103  0  0.00    651
#> 104  0  0.00    622
#> 105  0  0.00    594
#> 106  0  0.00    567
#> 107  0  0.00    541
#> 108  0  0.00    517
#> 109  0  0.00    493
#> 110  0  0.00    471
#> 111  0  0.00    450
#> 112  0  0.00    429
#> 113  0  0.00    410
#> 114  0  0.00    391
#> 115  0  0.00    374
#> 116  0  0.00    357
#> 117  0  0.00    341
#> 118  0  0.00    325
#> 119  0  0.00    310
#> 120  0  0.00    296
#> 121  0  0.00    283
#> 122  0  0.00    270
#> 123  0  0.00    258
#> 124  0  0.00    246
#> 125  0  0.00    235
#> 126  0  0.00    225
#> 127  0  0.00    214
#> 128  0  0.00    205
#> 129  0  0.00    196
#> 130  0  0.00    187
#> 131  0  0.00    178
#> 132  0  0.00    170
#> 133  0  0.00    162
#> 134  0  0.00    155
#> 135  0  0.00    148
#> 136  0  0.00    141
#> 137  0  0.00    135
#> 138  0  0.00    129
#> 139  0  0.00    123
#> 140  1  0.10    118
#> 141  1  0.77    112
#> 142  1  1.39    107
#> 143  1  1.97    102
#> 144  1  2.52     98
#> 145  1  3.03     93
#> 146  1  3.51     89
#> 147  1  3.97     85
#> 148  1  4.39     81
#> 149  1  4.80     77
#> 150  1  5.18     74
#> 151  1  5.54     71
#> 152  1  5.87     67
#> 153  1  6.20     64
#> 154  1  6.50     61
#> 155  1  6.79     59
#> 156  1  7.06     56
#> 157  1  7.32     53
#> 158  1  7.56     51
#> 159  1  7.79     49
#> 160  1  8.01     47
#> 161  1  8.21     44
#> 162  1  8.41     42
#> 163  1  8.59     41
#> 164  1  8.77     39
#> 165  1  8.93     37
#> 166  1  9.09     35
#> 167  1  9.24     34
#> 168  1  9.38     32
#> 169  1  9.51     31
#> 170  1  9.64     29
#> 171  1  9.75     28
#> 172  1  9.86     27
#> 173  1  9.97     26
#> 174  1 10.06     24
#> 175  1 10.16     23
#> 176  1 10.24     22
#> 177  1 10.33     21
#> 178  1 10.40     20
#> 179  1 10.47     19
#> 180  1 10.54     18
#> 181  1 10.61     18
#> 182  1 10.66     17
#> 183  1 10.72     16
#> 184  1 10.77     15
#> 185  1 10.82     15
#> 186  1 10.86     14
#> 187  1 10.91     13
#> 188  1 10.95     13
#> 189  2 11.05     12
#> 190  2 11.16     12
#> 191  2 11.25     11
#> 192  2 11.34     11
#> 193  2 11.42     10
#> 194  2 11.50     10
#> 195  2 11.57      9
#> 196  2 11.63      9
#> 197  2 11.69      8
#> 198  2 11.74      8
#> 199  2 11.79      8
#> 200  2 11.83      7
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev  Lambda
#> 1    0  0.00 26480.0
#> 2    0  0.00 25280.0
#> 3    0  0.00 24140.0
#> 4    0  0.00 23050.0
#> 5    0  0.00 22000.0
#> 6    0  0.00 21010.0
#> 7    0  0.00 20060.0
#> 8    0  0.00 19150.0
#> 9    0  0.00 18280.0
#> 10   0  0.00 17460.0
#> 11   0  0.00 16670.0
#> 12   0  0.00 15910.0
#> 13   0  0.00 15190.0
#> 14   0  0.00 14510.0
#> 15   0  0.00 13850.0
#> 16   0  0.00 13220.0
#> 17   0  0.00 12630.0
#> 18   0  0.00 12060.0
#> 19   0  0.00 11510.0
#> 20   0  0.00 10990.0
#> 21   0  0.00 10490.0
#> 22   0  0.00 10020.0
#> 23   0  0.00  9565.0
#> 24   0  0.00  9132.0
#> 25   0  0.00  8719.0
#> 26   0  0.00  8325.0
#> 27   0  0.00  7948.0
#> 28   0  0.00  7589.0
#> 29   0  0.00  7246.0
#> 30   0  0.00  6918.0
#> 31   0  0.00  6605.0
#> 32   0  0.00  6306.0
#> 33   0  0.00  6021.0
#> 34   0  0.00  5749.0
#> 35   0  0.00  5489.0
#> 36   0  0.00  5240.0
#> 37   0  0.00  5003.0
#> 38   0  0.00  4777.0
#> 39   0  0.00  4561.0
#> 40   0  0.00  4355.0
#> 41   0  0.00  4158.0
#> 42   0  0.00  3970.0
#> 43   0  0.00  3790.0
#> 44   0  0.00  3619.0
#> 45   0  0.00  3455.0
#> 46   0  0.00  3299.0
#> 47   0  0.00  3150.0
#> 48   0  0.00  3007.0
#> 49   0  0.00  2871.0
#> 50   0  0.00  2741.0
#> 51   0  0.00  2617.0
#> 52   0  0.00  2499.0
#> 53   0  0.00  2386.0
#> 54   0  0.00  2278.0
#> 55   0  0.00  2175.0
#> 56   0  0.00  2077.0
#> 57   0  0.00  1983.0
#> 58   0  0.00  1893.0
#> 59   0  0.00  1807.0
#> 60   0  0.00  1726.0
#> 61   0  0.00  1648.0
#> 62   0  0.00  1573.0
#> 63   0  0.00  1502.0
#> 64   0  0.00  1434.0
#> 65   0  0.00  1369.0
#> 66   0  0.00  1307.0
#> 67   0  0.00  1248.0
#> 68   0  0.00  1192.0
#> 69   0  0.00  1138.0
#> 70   0  0.00  1086.0
#> 71   0  0.00  1037.0
#> 72   0  0.00   990.3
#> 73   0  0.00   945.5
#> 74   0  0.00   902.7
#> 75   0  0.00   861.9
#> 76   0  0.00   822.9
#> 77   0  0.00   785.7
#> 78   0  0.00   750.2
#> 79   0  0.00   716.2
#> 80   0  0.00   683.8
#> 81   0  0.00   652.9
#> 82   0  0.00   623.4
#> 83   0  0.00   595.2
#> 84   0  0.00   568.3
#> 85   0  0.00   542.6
#> 86   0  0.00   518.0
#> 87   0  0.00   494.6
#> 88   0  0.00   472.2
#> 89   0  0.00   450.9
#> 90   0  0.00   430.5
#> 91   0  0.00   411.0
#> 92   0  0.00   392.4
#> 93   0  0.00   374.7
#> 94   0  0.00   357.7
#> 95   0  0.00   341.5
#> 96   0  0.00   326.1
#> 97   0  0.00   311.3
#> 98   0  0.00   297.3
#> 99   0  0.00   283.8
#> 100  0  0.00   271.0
#> 101  0  0.00   258.7
#> 102  0  0.00   247.0
#> 103  0  0.00   235.9
#> 104  0  0.00   225.2
#> 105  1  0.55   215.0
#> 106  1  1.35   205.3
#> 107  1  2.09   196.0
#> 108  1  2.78   187.1
#> 109  1  3.42   178.7
#> 110  1  4.02   170.6
#> 111  1  4.59   162.9
#> 112  1  5.12   155.5
#> 113  1  5.62   148.5
#> 114  1  6.10   141.8
#> 115  1  6.54   135.3
#> 116  1  6.96   129.2
#> 117  1  7.36   123.4
#> 118  1  7.73   117.8
#> 119  1  8.09   112.5
#> 120  1  8.43   107.4
#> 121  1  8.76   102.5
#> 122  1  9.07    97.9
#> 123  1  9.36    93.5
#> 124  1  9.64    89.2
#> 125  1  9.91    85.2
#> 126  1 10.17    81.3
#> 127  1 10.41    77.7
#> 128  1 10.65    74.2
#> 129  1 10.87    70.8
#> 130  1 11.08    67.6
#> 131  1 11.29    64.5
#> 132  1 11.49    61.6
#> 133  1 11.68    58.8
#> 134  1 11.86    56.2
#> 135  1 12.03    53.6
#> 136  1 12.20    51.2
#> 137  1 12.36    48.9
#> 138  1 12.51    46.7
#> 139  1 12.66    44.6
#> 140  1 12.80    42.5
#> 141  1 12.94    40.6
#> 142  1 13.07    38.8
#> 143  1 13.20    37.0
#> 144  1 13.32    35.4
#> 145  1 13.44    33.8
#> 146  1 13.55    32.2
#> 147  1 13.66    30.8
#> 148  1 13.76    29.4
#> 149  1 13.86    28.1
#> 150  1 13.96    26.8
#> 151  1 14.05    25.6
#> 152  1 14.14    24.4
#> 153  1 14.22    23.3
#> 154  1 14.31    22.3
#> 155  1 14.38    21.2
#> 156  1 14.46    20.3
#> 157  1 14.53    19.4
#> 158  1 14.60    18.5
#> 159  1 14.66    17.7
#> 160  1 14.73    16.9
#> 161  1 14.79    16.1
#> 162  1 14.85    15.4
#> 163  1 14.90    14.7
#> 164  1 14.95    14.0
#> 165  1 15.00    13.4
#> 166  1 15.05    12.8
#> 167  1 15.09    12.2
#> 168  1 15.14    11.6
#> 169  1 15.18    11.1
#> 170  1 15.22    10.6
#> 171  1 15.25    10.1
#> 172  1 15.29     9.7
#> 173  1 15.32     9.2
#> 174  1 15.35     8.8
#> 175  1 15.38     8.4
#> 176  2 15.48     8.0
#> 177  2 15.82     7.7
#> 178  2 16.13     7.3
#> 179  2 16.41     7.0
#> 180  2 16.67     6.7
#> 181  2 16.90     6.4
#> 182  2 17.11     6.1
#> 183  2 17.30     5.8
#> 184  2 17.47     5.6
#> 185  2 17.63     5.3
#> 186  2 17.78     5.1
#> 187  2 17.91     4.8
#> 188  2 18.03     4.6
#> 189  2 18.14     4.4
#> 190  2 18.24     4.2
#> 191  2 18.33     4.0
#> 192  2 18.41     3.8
#> 193  2 18.49     3.7
#> 194  2 18.56     3.5
#> 195  2 18.63     3.3
#> 196  2 18.69     3.2
#> 197  2 18.74     3.0
#> 198  2 18.79     2.9
#> 199  2 18.83     2.8
#> 200  2 18.88     2.6
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev  Lambda
#> 1    0  0.00 117.600
#> 2    0  0.00 112.200
#> 3    0  0.00 107.200
#> 4    0  0.00 102.300
#> 5    0  0.00  97.690
#> 6    0  0.00  93.270
#> 7    0  0.00  89.050
#> 8    0  0.00  85.020
#> 9    0  0.00  81.180
#> 10   0  0.00  77.510
#> 11   0  0.00  74.000
#> 12   0  0.00  70.650
#> 13   0  0.00  67.460
#> 14   0  0.00  64.410
#> 15   0  0.00  61.490
#> 16   0  0.00  58.710
#> 17   0  0.00  56.060
#> 18   0  0.00  53.520
#> 19   0  0.00  51.100
#> 20   0  0.00  48.790
#> 21   0  0.00  46.580
#> 22   0  0.00  44.480
#> 23   0  0.00  42.470
#> 24   0  0.00  40.540
#> 25   0  0.00  38.710
#> 26   0  0.00  36.960
#> 27   0  0.00  35.290
#> 28   0  0.00  33.690
#> 29   0  0.00  32.170
#> 30   0  0.00  30.710
#> 31   0  0.00  29.320
#> 32   0  0.00  28.000
#> 33   0  0.00  26.730
#> 34   0  0.00  25.520
#> 35   0  0.00  24.370
#> 36   0  0.00  23.270
#> 37   0  0.00  22.210
#> 38   0  0.00  21.210
#> 39   0  0.00  20.250
#> 40   0  0.00  19.330
#> 41   0  0.00  18.460
#> 42   0  0.00  17.620
#> 43   0  0.00  16.830
#> 44   0  0.00  16.070
#> 45   0  0.00  15.340
#> 46   0  0.00  14.650
#> 47   0  0.00  13.980
#> 48   0  0.00  13.350
#> 49   0  0.00  12.750
#> 50   0  0.00  12.170
#> 51   0  0.00  11.620
#> 52   0  0.00  11.090
#> 53   0  0.00  10.590
#> 54   0  0.00  10.110
#> 55   0  0.00   9.656
#> 56   0  0.00   9.220
#> 57   0  0.00   8.803
#> 58   0  0.00   8.405
#> 59   0  0.00   8.024
#> 60   0  0.00   7.662
#> 61   0  0.00   7.315
#> 62   0  0.00   6.984
#> 63   0  0.00   6.668
#> 64   0  0.00   6.367
#> 65   0  0.00   6.079
#> 66   0  0.00   5.804
#> 67   0  0.00   5.541
#> 68   0  0.00   5.291
#> 69   0  0.00   5.051
#> 70   0  0.00   4.823
#> 71   0  0.00   4.605
#> 72   0  0.00   4.397
#> 73   0  0.00   4.198
#> 74   0  0.00   4.008
#> 75   0  0.00   3.827
#> 76   0  0.00   3.653
#> 77   0  0.00   3.488
#> 78   0  0.00   3.330
#> 79   0  0.00   3.180
#> 80   0  0.00   3.036
#> 81   0  0.00   2.899
#> 82   0  0.00   2.768
#> 83   0  0.00   2.642
#> 84   0  0.00   2.523
#> 85   0  0.00   2.409
#> 86   0  0.00   2.300
#> 87   0  0.00   2.196
#> 88   0  0.00   2.097
#> 89   0  0.00   2.002
#> 90   0  0.00   1.911
#> 91   0  0.00   1.825
#> 92   0  0.00   1.742
#> 93   0  0.00   1.663
#> 94   0  0.00   1.588
#> 95   0  0.00   1.516
#> 96   0  0.00   1.448
#> 97   0  0.00   1.382
#> 98   0  0.00   1.320
#> 99   0  0.00   1.260
#> 100  0  0.00   1.203
#> 101  0  0.00   1.149
#> 102  0  0.00   1.097
#> 103  0  0.00   1.047
#> 104  0  0.00   1.000
#> 105  1  0.55   0.955
#> 106  1  1.35   0.911
#> 107  1  2.09   0.870
#> 108  1  2.78   0.831
#> 109  1  3.42   0.793
#> 110  1  4.02   0.757
#> 111  1  4.59   0.723
#> 112  1  5.12   0.690
#> 113  1  5.62   0.659
#> 114  1  6.10   0.629
#> 115  1  6.54   0.601
#> 116  1  6.96   0.574
#> 117  1  7.36   0.548
#> 118  1  7.73   0.523
#> 119  1  8.09   0.499
#> 120  1  8.43   0.477
#> 121  1  8.76   0.455
#> 122  1  9.07   0.435
#> 123  1  9.36   0.415
#> 124  1  9.64   0.396
#> 125  1  9.91   0.378
#> 126  1 10.17   0.361
#> 127  1 10.41   0.345
#> 128  1 10.65   0.329
#> 129  1 10.87   0.314
#> 130  1 11.08   0.300
#> 131  1 11.29   0.286
#> 132  1 11.49   0.274
#> 133  1 11.68   0.261
#> 134  1 11.86   0.249
#> 135  1 12.03   0.238
#> 136  1 12.20   0.227
#> 137  1 12.36   0.217
#> 138  1 12.51   0.207
#> 139  1 12.66   0.198
#> 140  1 12.80   0.189
#> 141  1 12.94   0.180
#> 142  1 13.07   0.172
#> 143  1 13.20   0.164
#> 144  1 13.32   0.157
#> 145  1 13.44   0.150
#> 146  1 13.55   0.143
#> 147  1 13.66   0.137
#> 148  1 13.76   0.130
#> 149  1 13.86   0.125
#> 150  1 13.96   0.119
#> 151  1 14.05   0.114
#> 152  1 14.14   0.108
#> 153  1 14.22   0.104
#> 154  1 14.31   0.099
#> 155  1 14.38   0.094
#> 156  1 14.46   0.090
#> 157  1 14.53   0.086
#> 158  1 14.60   0.082
#> 159  1 14.66   0.078
#> 160  1 14.73   0.075
#> 161  1 14.79   0.071
#> 162  1 14.85   0.068
#> 163  1 14.90   0.065
#> 164  1 14.95   0.062
#> 165  1 15.00   0.059
#> 166  1 15.05   0.057
#> 167  1 15.09   0.054
#> 168  1 15.14   0.052
#> 169  1 15.18   0.049
#> 170  1 15.22   0.047
#> 171  1 15.25   0.045
#> 172  1 15.29   0.043
#> 173  1 15.32   0.041
#> 174  1 15.35   0.039
#> 175  1 15.38   0.037
#> 176  1 15.41   0.036
#> 177  1 15.43   0.034
#> 178  1 15.46   0.033
#> 179  1 15.48   0.031
#> 180  1 15.50   0.030
#> 181  1 15.52   0.028
#> 182  1 15.54   0.027
#> 183  1 15.56   0.026
#> 184  1 15.58   0.025
#> 185  1 15.59   0.024
#> 186  1 15.61   0.022
#> 187  1 15.62   0.021
#> 188  1 15.63   0.020
#> 189  1 15.65   0.020
#> 190  1 15.66   0.019
#> 191  1 15.67   0.018
#> 192  1 15.68   0.017
#> 193  1 15.69   0.016
#> 194  1 15.70   0.016
#> 195  1 15.70   0.015
#> 196  1 15.71   0.014
#> 197  1 15.72   0.014
#> 198  1 15.73   0.013
#> 199  1 15.73   0.012
#> 200  1 15.74   0.012
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev  Lambda
#> 1    0  0.00 11.2200
#> 2    0  0.00 10.7100
#> 3    0  0.00 10.2300
#> 4    0  0.00  9.7630
#> 5    0  0.00  9.3210
#> 6    0  0.00  8.8990
#> 7    0  0.00  8.4970
#> 8    0  0.00  8.1130
#> 9    0  0.00  7.7460
#> 10   0  0.00  7.3950
#> 11   0  0.00  7.0610
#> 12   0  0.00  6.7420
#> 13   0  0.00  6.4370
#> 14   0  0.00  6.1460
#> 15   0  0.00  5.8680
#> 16   0  0.00  5.6020
#> 17   0  0.00  5.3490
#> 18   0  0.00  5.1070
#> 19   0  0.00  4.8760
#> 20   0  0.00  4.6550
#> 21   0  0.00  4.4450
#> 22   0  0.00  4.2440
#> 23   0  0.00  4.0520
#> 24   0  0.00  3.8690
#> 25   0  0.00  3.6940
#> 26   0  0.00  3.5270
#> 27   0  0.00  3.3670
#> 28   0  0.00  3.2150
#> 29   0  0.00  3.0690
#> 30   0  0.00  2.9310
#> 31   0  0.00  2.7980
#> 32   0  0.00  2.6710
#> 33   0  0.00  2.5510
#> 34   0  0.00  2.4350
#> 35   0  0.00  2.3250
#> 36   0  0.00  2.2200
#> 37   0  0.00  2.1200
#> 38   0  0.00  2.0240
#> 39   0  0.00  1.9320
#> 40   0  0.00  1.8450
#> 41   0  0.00  1.7610
#> 42   0  0.00  1.6820
#> 43   0  0.00  1.6060
#> 44   0  0.00  1.5330
#> 45   0  0.00  1.4640
#> 46   0  0.00  1.3970
#> 47   0  0.00  1.3340
#> 48   0  0.00  1.2740
#> 49   0  0.00  1.2160
#> 50   0  0.00  1.1610
#> 51   0  0.00  1.1090
#> 52   0  0.00  1.0590
#> 53   0  0.00  1.0110
#> 54   0  0.00  0.9650
#> 55   0  0.00  0.9214
#> 56   0  0.00  0.8797
#> 57   0  0.00  0.8399
#> 58   0  0.00  0.8019
#> 59   0  0.00  0.7657
#> 60   0  0.00  0.7310
#> 61   0  0.00  0.6980
#> 62   0  0.00  0.6664
#> 63   0  0.00  0.6363
#> 64   0  0.00  0.6075
#> 65   0  0.00  0.5800
#> 66   0  0.00  0.5538
#> 67   0  0.00  0.5287
#> 68   0  0.00  0.5048
#> 69   0  0.00  0.4820
#> 70   0  0.00  0.4602
#> 71   0  0.00  0.4394
#> 72   0  0.00  0.4195
#> 73   0  0.00  0.4005
#> 74   0  0.00  0.3824
#> 75   0  0.00  0.3651
#> 76   0  0.00  0.3486
#> 77   0  0.00  0.3328
#> 78   0  0.00  0.3178
#> 79   0  0.00  0.3034
#> 80   0  0.00  0.2897
#> 81   0  0.00  0.2766
#> 82   0  0.00  0.2641
#> 83   0  0.00  0.2521
#> 84   0  0.00  0.2407
#> 85   0  0.00  0.2298
#> 86   0  0.00  0.2194
#> 87   0  0.00  0.2095
#> 88   0  0.00  0.2000
#> 89   0  0.00  0.1910
#> 90   0  0.00  0.1824
#> 91   0  0.00  0.1741
#> 92   0  0.00  0.1662
#> 93   0  0.00  0.1587
#> 94   0  0.00  0.1515
#> 95   0  0.00  0.1447
#> 96   0  0.00  0.1381
#> 97   0  0.00  0.1319
#> 98   0  0.00  0.1259
#> 99   0  0.00  0.1202
#> 100  0  0.00  0.1148
#> 101  0  0.00  0.1096
#> 102  0  0.00  0.1046
#> 103  0  0.00  0.0999
#> 104  0  0.00  0.0954
#> 105  1  0.55  0.0911
#> 106  1  1.35  0.0870
#> 107  1  2.09  0.0830
#> 108  1  2.78  0.0793
#> 109  1  3.42  0.0757
#> 110  1  4.02  0.0723
#> 111  1  4.59  0.0690
#> 112  1  5.12  0.0659
#> 113  1  5.62  0.0629
#> 114  1  6.10  0.0600
#> 115  1  6.54  0.0573
#> 116  1  6.96  0.0547
#> 117  1  7.36  0.0523
#> 118  1  7.73  0.0499
#> 119  1  8.09  0.0476
#> 120  1  8.43  0.0455
#> 121  1  8.76  0.0434
#> 122  1  9.07  0.0415
#> 123  1  9.36  0.0396
#> 124  1  9.64  0.0378
#> 125  1  9.91  0.0361
#> 126  1 10.17  0.0345
#> 127  1 10.41  0.0329
#> 128  1 10.65  0.0314
#> 129  1 10.87  0.0300
#> 130  1 11.08  0.0286
#> 131  1 11.29  0.0273
#> 132  1 11.49  0.0261
#> 133  1 11.68  0.0249
#> 134  1 11.86  0.0238
#> 135  1 12.03  0.0227
#> 136  1 12.20  0.0217
#> 137  1 12.36  0.0207
#> 138  1 12.51  0.0198
#> 139  1 12.66  0.0189
#> 140  1 12.80  0.0180
#> 141  1 12.94  0.0172
#> 142  1 13.07  0.0164
#> 143  1 13.20  0.0157
#> 144  1 13.32  0.0150
#> 145  1 13.44  0.0143
#> 146  1 13.55  0.0137
#> 147  1 13.66  0.0130
#> 148  1 13.76  0.0124
#> 149  1 13.86  0.0119
#> 150  1 13.96  0.0114
#> 151  1 14.05  0.0108
#> 152  1 14.14  0.0103
#> 153  1 14.22  0.0099
#> 154  1 14.31  0.0094
#> 155  1 14.38  0.0090
#> 156  1 14.46  0.0086
#> 157  1 14.53  0.0082
#> 158  1 14.60  0.0078
#> 159  1 14.66  0.0075
#> 160  1 14.73  0.0071
#> 161  1 14.79  0.0068
#> 162  1 14.85  0.0065
#> 163  1 14.90  0.0062
#> 164  1 14.95  0.0059
#> 165  1 15.00  0.0057
#> 166  1 15.05  0.0054
#> 167  1 15.09  0.0052
#> 168  1 15.14  0.0049
#> 169  1 15.18  0.0047
#> 170  1 15.22  0.0045
#> 171  1 15.25  0.0043
#> 172  1 15.29  0.0041
#> 173  1 15.32  0.0039
#> 174  1 15.35  0.0037
#> 175  1 15.38  0.0036
#> 176  1 15.41  0.0034
#> 177  1 15.43  0.0033
#> 178  1 15.46  0.0031
#> 179  1 15.48  0.0030
#> 180  1 15.50  0.0028
#> 181  1 15.52  0.0027
#> 182  1 15.54  0.0026
#> 183  1 15.56  0.0025
#> 184  1 15.58  0.0024
#> 185  1 15.59  0.0022
#> 186  1 15.61  0.0021
#> 187  1 15.62  0.0020
#> 188  1 15.63  0.0020
#> 189  1 15.65  0.0019
#> 190  1 15.66  0.0018
#> 191  1 15.67  0.0017
#> 192  1 15.68  0.0016
#> 193  1 15.69  0.0016
#> 194  1 15.70  0.0015
#> 195  1 15.70  0.0014
#> 196  1 15.71  0.0014
#> 197  1 15.72  0.0013
#> 198  1 15.73  0.0012
#> 199  1 15.73  0.0012
#> 200  1 15.74  0.0011
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev  Lambda
#> 1    0  0.00 2.41100
#> 2    0  0.00 2.30200
#> 3    0  0.00 2.19800
#> 4    0  0.00 2.09900
#> 5    0  0.00 2.00400
#> 6    0  0.00 1.91300
#> 7    0  0.00 1.82700
#> 8    0  0.00 1.74400
#> 9    0  0.00 1.66500
#> 10   0  0.00 1.59000
#> 11   0  0.00 1.51800
#> 12   0  0.00 1.44900
#> 13   0  0.00 1.38400
#> 14   0  0.00 1.32100
#> 15   0  0.00 1.26100
#> 16   0  0.00 1.20400
#> 17   0  0.00 1.15000
#> 18   0  0.00 1.09800
#> 19   0  0.00 1.04800
#> 20   0  0.00 1.00100
#> 21   0  0.00 0.95550
#> 22   0  0.00 0.91230
#> 23   0  0.00 0.87100
#> 24   0  0.00 0.83160
#> 25   0  0.00 0.79400
#> 26   0  0.00 0.75810
#> 27   0  0.00 0.72380
#> 28   0  0.00 0.69110
#> 29   0  0.00 0.65980
#> 30   0  0.00 0.63000
#> 31   0  0.00 0.60150
#> 32   0  0.00 0.57430
#> 33   0  0.00 0.54830
#> 34   0  0.00 0.52350
#> 35   0  0.00 0.49980
#> 36   0  0.00 0.47720
#> 37   0  0.00 0.45560
#> 38   0  0.00 0.43500
#> 39   0  0.00 0.41530
#> 40   0  0.00 0.39660
#> 41   0  0.00 0.37860
#> 42   0  0.00 0.36150
#> 43   0  0.00 0.34510
#> 44   0  0.00 0.32950
#> 45   0  0.00 0.31460
#> 46   0  0.00 0.30040
#> 47   0  0.00 0.28680
#> 48   0  0.00 0.27380
#> 49   0  0.00 0.26150
#> 50   0  0.00 0.24960
#> 51   0  0.00 0.23830
#> 52   0  0.00 0.22760
#> 53   0  0.00 0.21730
#> 54   0  0.00 0.20740
#> 55   0  0.00 0.19810
#> 56   0  0.00 0.18910
#> 57   0  0.00 0.18060
#> 58   0  0.00 0.17240
#> 59   0  0.00 0.16460
#> 60   0  0.00 0.15710
#> 61   0  0.00 0.15000
#> 62   0  0.00 0.14330
#> 63   0  0.00 0.13680
#> 64   0  0.00 0.13060
#> 65   0  0.00 0.12470
#> 66   0  0.00 0.11900
#> 67   0  0.00 0.11370
#> 68   0  0.00 0.10850
#> 69   0  0.00 0.10360
#> 70   0  0.00 0.09892
#> 71   0  0.00 0.09445
#> 72   0  0.00 0.09018
#> 73   0  0.00 0.08610
#> 74   0  0.00 0.08220
#> 75   0  0.00 0.07849
#> 76   0  0.00 0.07494
#> 77   0  0.00 0.07155
#> 78   0  0.00 0.06831
#> 79   0  0.00 0.06522
#> 80   0  0.00 0.06227
#> 81   0  0.00 0.05946
#> 82   0  0.00 0.05677
#> 83   0  0.00 0.05420
#> 84   0  0.00 0.05175
#> 85   0  0.00 0.04941
#> 86   0  0.00 0.04717
#> 87   0  0.00 0.04504
#> 88   0  0.00 0.04300
#> 89   0  0.00 0.04106
#> 90   0  0.00 0.03920
#> 91   0  0.00 0.03743
#> 92   0  0.00 0.03573
#> 93   0  0.00 0.03412
#> 94   0  0.00 0.03257
#> 95   0  0.00 0.03110
#> 96   0  0.00 0.02969
#> 97   0  0.00 0.02835
#> 98   0  0.00 0.02707
#> 99   0  0.00 0.02585
#> 100  0  0.00 0.02468
#> 101  0  0.00 0.02356
#> 102  0  0.00 0.02249
#> 103  0  0.00 0.02148
#> 104  0  0.00 0.02051
#> 105  1  0.55 0.01958
#> 106  1  1.35 0.01869
#> 107  1  2.09 0.01785
#> 108  1  2.78 0.01704
#> 109  1  3.42 0.01627
#> 110  1  4.02 0.01553
#> 111  1  4.59 0.01483
#> 112  1  5.12 0.01416
#> 113  1  5.62 0.01352
#> 114  1  6.10 0.01291
#> 115  1  6.54 0.01232
#> 116  1  6.96 0.01177
#> 117  1  7.36 0.01123
#> 118  1  7.73 0.01073
#> 119  1  8.09 0.01024
#> 120  1  8.43 0.00978
#> 121  1  8.76 0.00934
#> 122  1  9.07 0.00891
#> 123  1  9.36 0.00851
#> 124  1  9.64 0.00813
#> 125  1  9.91 0.00776
#> 126  1 10.17 0.00741
#> 127  1 10.41 0.00707
#> 128  1 10.65 0.00675
#> 129  1 10.87 0.00645
#> 130  1 11.08 0.00616
#> 131  1 11.29 0.00588
#> 132  1 11.49 0.00561
#> 133  1 11.68 0.00536
#> 134  1 11.86 0.00511
#> 135  1 12.03 0.00488
#> 136  1 12.20 0.00466
#> 137  1 12.36 0.00445
#> 138  1 12.51 0.00425
#> 139  1 12.66 0.00406
#> 140  1 12.80 0.00387
#> 141  1 12.94 0.00370
#> 142  1 13.07 0.00353
#> 143  1 13.20 0.00337
#> 144  1 13.32 0.00322
#> 145  1 13.44 0.00307
#> 146  1 13.55 0.00294
#> 147  1 13.66 0.00280
#> 148  1 13.76 0.00268
#> 149  1 13.86 0.00255
#> 150  1 13.96 0.00244
#> 151  1 14.05 0.00233
#> 152  1 14.14 0.00222
#> 153  1 14.22 0.00212
#> 154  1 14.31 0.00203
#> 155  1 14.38 0.00194
#> 156  1 14.46 0.00185
#> 157  1 14.53 0.00176
#> 158  1 14.60 0.00168
#> 159  1 14.66 0.00161
#> 160  1 14.73 0.00154
#> 161  1 14.79 0.00147
#> 162  1 14.85 0.00140
#> 163  1 14.90 0.00134
#> 164  1 14.95 0.00128
#> 165  1 15.00 0.00122
#> 166  1 15.05 0.00116
#> 167  1 15.09 0.00111
#> 168  1 15.14 0.00106
#> 169  1 15.18 0.00101
#> 170  1 15.22 0.00097
#> 171  1 15.25 0.00092
#> 172  1 15.29 0.00088
#> 173  1 15.32 0.00084
#> 174  1 15.35 0.00080
#> 175  1 15.38 0.00077
#> 176  1 15.41 0.00073
#> 177  1 15.43 0.00070
#> 178  1 15.46 0.00067
#> 179  1 15.48 0.00064
#> 180  1 15.50 0.00061
#> 181  1 15.52 0.00058
#> 182  1 15.54 0.00055
#> 183  1 15.56 0.00053
#> 184  1 15.58 0.00051
#> 185  1 15.59 0.00048
#> 186  1 15.61 0.00046
#> 187  1 15.62 0.00044
#> 188  1 15.63 0.00042
#> 189  1 15.65 0.00040
#> 190  1 15.66 0.00038
#> 191  1 15.67 0.00037
#> 192  1 15.68 0.00035
#> 193  1 15.69 0.00033
#> 194  1 15.70 0.00032
#> 195  1 15.70 0.00030
#> 196  1 15.71 0.00029
#> 197  1 15.72 0.00028
#> 198  1 15.73 0.00026
#> 199  1 15.73 0.00025
#> 200  1 15.74 0.00024
#> 
#> $`1`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev Lambda
#> 1    0  0.00 8699.0
#> 2    0  0.00 8306.0
#> 3    0  0.00 7930.0
#> 4    0  0.00 7571.0
#> 5    0  0.00 7229.0
#> 6    0  0.00 6902.0
#> 7    0  0.00 6590.0
#> 8    0  0.00 6292.0
#> 9    0  0.00 6007.0
#> 10   0  0.00 5735.0
#> 11   0  0.00 5476.0
#> 12   0  0.00 5228.0
#> 13   0  0.00 4992.0
#> 14   0  0.00 4766.0
#> 15   0  0.00 4551.0
#> 16   0  0.00 4345.0
#> 17   0  0.00 4148.0
#> 18   0  0.00 3961.0
#> 19   0  0.00 3781.0
#> 20   0  0.00 3610.0
#> 21   0  0.00 3447.0
#> 22   0  0.00 3291.0
#> 23   0  0.00 3142.0
#> 24   0  0.00 3000.0
#> 25   0  0.00 2865.0
#> 26   0  0.00 2735.0
#> 27   0  0.00 2611.0
#> 28   0  0.00 2493.0
#> 29   0  0.00 2380.0
#> 30   0  0.00 2273.0
#> 31   0  0.00 2170.0
#> 32   0  0.00 2072.0
#> 33   0  0.00 1978.0
#> 34   0  0.00 1889.0
#> 35   0  0.00 1803.0
#> 36   0  0.00 1722.0
#> 37   0  0.00 1644.0
#> 38   0  0.00 1569.0
#> 39   0  0.00 1498.0
#> 40   0  0.00 1431.0
#> 41   0  0.00 1366.0
#> 42   0  0.00 1304.0
#> 43   0  0.00 1245.0
#> 44   0  0.00 1189.0
#> 45   0  0.00 1135.0
#> 46   0  0.00 1084.0
#> 47   0  0.00 1035.0
#> 48   0  0.00  988.0
#> 49   0  0.00  943.3
#> 50   0  0.00  900.6
#> 51   0  0.00  859.9
#> 52   0  0.00  821.0
#> 53   0  0.00  783.9
#> 54   0  0.00  748.4
#> 55   0  0.00  714.6
#> 56   0  0.00  682.2
#> 57   0  0.00  651.4
#> 58   0  0.00  621.9
#> 59   0  0.00  593.8
#> 60   0  0.00  566.9
#> 61   0  0.00  541.3
#> 62   0  0.00  516.8
#> 63   0  0.00  493.4
#> 64   0  0.00  471.1
#> 65   0  0.00  449.8
#> 66   0  0.00  429.5
#> 67   0  0.00  410.1
#> 68   0  0.00  391.5
#> 69   0  0.00  373.8
#> 70   0  0.00  356.9
#> 71   0  0.00  340.7
#> 72   0  0.00  325.3
#> 73   0  0.00  310.6
#> 74   0  0.00  296.6
#> 75   0  0.00  283.2
#> 76   0  0.00  270.4
#> 77   0  0.00  258.1
#> 78   0  0.00  246.5
#> 79   0  0.00  235.3
#> 80   0  0.00  224.7
#> 81   0  0.00  214.5
#> 82   0  0.00  204.8
#> 83   0  0.00  195.5
#> 84   0  0.00  186.7
#> 85   0  0.00  178.2
#> 86   0  0.00  170.2
#> 87   0  0.00  162.5
#> 88   0  0.00  155.1
#> 89   0  0.00  148.1
#> 90   0  0.00  141.4
#> 91   0  0.00  135.0
#> 92   0  0.00  128.9
#> 93   0  0.00  123.1
#> 94   0  0.00  117.5
#> 95   0  0.00  112.2
#> 96   0  0.00  107.1
#> 97   0  0.00  102.3
#> 98   0  0.00   97.7
#> 99   0  0.00   93.2
#> 100  0  0.00   89.0
#> 101  0  0.00   85.0
#> 102  0  0.00   81.2
#> 103  0  0.00   77.5
#> 104  0  0.00   74.0
#> 105  1  0.55   70.6
#> 106  1  1.35   67.4
#> 107  1  2.09   64.4
#> 108  1  2.78   61.5
#> 109  1  3.42   58.7
#> 110  1  4.02   56.0
#> 111  1  4.59   53.5
#> 112  1  5.12   51.1
#> 113  1  5.62   48.8
#> 114  1  6.10   46.6
#> 115  1  6.54   44.5
#> 116  1  6.96   42.5
#> 117  1  7.36   40.5
#> 118  1  7.73   38.7
#> 119  1  8.09   37.0
#> 120  1  8.43   35.3
#> 121  1  8.76   33.7
#> 122  1  9.07   32.2
#> 123  1  9.36   30.7
#> 124  1  9.64   29.3
#> 125  1  9.91   28.0
#> 126  1 10.17   26.7
#> 127  1 10.41   25.5
#> 128  1 10.65   24.4
#> 129  1 10.87   23.3
#> 130  1 11.08   22.2
#> 131  1 11.29   21.2
#> 132  1 11.49   20.2
#> 133  1 11.68   19.3
#> 134  1 11.86   18.4
#> 135  1 12.03   17.6
#> 136  1 12.20   16.8
#> 137  1 12.36   16.1
#> 138  1 12.51   15.3
#> 139  1 12.66   14.6
#> 140  1 12.80   14.0
#> 141  1 12.94   13.3
#> 142  1 13.07   12.7
#> 143  1 13.20   12.2
#> 144  1 13.32   11.6
#> 145  1 13.44   11.1
#> 146  1 13.55   10.6
#> 147  1 13.66   10.1
#> 148  1 13.76    9.7
#> 149  1 13.86    9.2
#> 150  1 13.96    8.8
#> 151  1 14.05    8.4
#> 152  1 14.14    8.0
#> 153  1 14.22    7.7
#> 154  1 14.31    7.3
#> 155  1 14.38    7.0
#> 156  1 14.46    6.7
#> 157  1 14.53    6.4
#> 158  1 14.60    6.1
#> 159  1 14.66    5.8
#> 160  1 14.73    5.5
#> 161  1 14.79    5.3
#> 162  1 14.85    5.0
#> 163  1 14.90    4.8
#> 164  1 14.95    4.6
#> 165  1 15.00    4.4
#> 166  1 15.05    4.2
#> 167  1 15.09    4.0
#> 168  1 15.14    3.8
#> 169  1 15.18    3.7
#> 170  1 15.22    3.5
#> 171  1 15.25    3.3
#> 172  1 15.29    3.2
#> 173  1 15.32    3.0
#> 174  1 15.35    2.9
#> 175  1 15.38    2.8
#> 176  1 15.41    2.6
#> 177  1 15.43    2.5
#> 178  1 15.46    2.4
#> 179  1 15.48    2.3
#> 180  1 15.50    2.2
#> 181  1 15.52    2.1
#> 182  2 15.69    2.0
#> 183  2 15.91    1.9
#> 184  2 16.10    1.8
#> 185  2 16.27    1.7
#> 186  2 16.43    1.7
#> 187  2 16.57    1.6
#> 188  2 16.70    1.5
#> 189  2 16.81    1.4
#> 190  2 16.92    1.4
#> 191  2 17.02    1.3
#> 192  2 17.10    1.3
#> 193  2 17.18    1.2
#> 194  2 17.26    1.1
#> 195  2 17.32    1.1
#> 196  2 17.38    1.0
#> 197  2 17.44    1.0
#> 198  2 17.49    1.0
#> 199  2 17.54    0.9
#> 200  2 17.58    0.9
#> 
#> $`0.230324074074074`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00 8.9330
#> 2    0 0.00 8.5290
#> 3    0 0.00 8.1440
#> 4    0 0.00 7.7750
#> 5    0 0.00 7.4240
#> 6    0 0.00 7.0880
#> 7    0 0.00 6.7670
#> 8    0 0.00 6.4610
#> 9    0 0.00 6.1690
#> 10   0 0.00 5.8900
#> 11   0 0.00 5.6240
#> 12   0 0.00 5.3690
#> 13   0 0.00 5.1260
#> 14   0 0.00 4.8950
#> 15   0 0.00 4.6730
#> 16   0 0.00 4.4620
#> 17   0 0.00 4.2600
#> 18   0 0.00 4.0670
#> 19   0 0.00 3.8830
#> 20   0 0.00 3.7080
#> 21   0 0.00 3.5400
#> 22   0 0.00 3.3800
#> 23   0 0.00 3.2270
#> 24   0 0.00 3.0810
#> 25   0 0.00 2.9420
#> 26   0 0.00 2.8090
#> 27   0 0.00 2.6820
#> 28   0 0.00 2.5600
#> 29   0 0.00 2.4450
#> 30   0 0.00 2.3340
#> 31   0 0.00 2.2280
#> 32   0 0.00 2.1280
#> 33   0 0.00 2.0310
#> 34   0 0.00 1.9400
#> 35   0 0.00 1.8520
#> 36   0 0.00 1.7680
#> 37   0 0.00 1.6880
#> 38   0 0.00 1.6120
#> 39   0 0.00 1.5390
#> 40   0 0.00 1.4690
#> 41   0 0.00 1.4030
#> 42   0 0.00 1.3390
#> 43   0 0.00 1.2790
#> 44   0 0.00 1.2210
#> 45   0 0.00 1.1660
#> 46   0 0.00 1.1130
#> 47   0 0.00 1.0630
#> 48   0 0.00 1.0150
#> 49   0 0.00 0.9687
#> 50   0 0.00 0.9249
#> 51   0 0.00 0.8831
#> 52   0 0.00 0.8431
#> 53   0 0.00 0.8050
#> 54   0 0.00 0.7686
#> 55   0 0.00 0.7338
#> 56   0 0.00 0.7006
#> 57   0 0.00 0.6689
#> 58   0 0.00 0.6387
#> 59   0 0.00 0.6098
#> 60   0 0.00 0.5822
#> 61   0 0.00 0.5559
#> 62   0 0.00 0.5307
#> 63   0 0.00 0.5067
#> 64   0 0.00 0.4838
#> 65   0 0.00 0.4619
#> 66   0 0.00 0.4410
#> 67   0 0.00 0.4211
#> 68   0 0.00 0.4021
#> 69   0 0.00 0.3839
#> 70   0 0.00 0.3665
#> 71   0 0.00 0.3499
#> 72   0 0.00 0.3341
#> 73   0 0.00 0.3190
#> 74   0 0.00 0.3046
#> 75   0 0.00 0.2908
#> 76   0 0.00 0.2776
#> 77   0 0.00 0.2651
#> 78   0 0.00 0.2531
#> 79   0 0.00 0.2416
#> 80   0 0.00 0.2307
#> 81   0 0.00 0.2203
#> 82   0 0.00 0.2103
#> 83   0 0.00 0.2008
#> 84   0 0.00 0.1917
#> 85   0 0.00 0.1831
#> 86   0 0.00 0.1748
#> 87   0 0.00 0.1669
#> 88   0 0.00 0.1593
#> 89   0 0.00 0.1521
#> 90   0 0.00 0.1452
#> 91   0 0.00 0.1387
#> 92   0 0.00 0.1324
#> 93   0 0.00 0.1264
#> 94   0 0.00 0.1207
#> 95   0 0.00 0.1152
#> 96   0 0.00 0.1100
#> 97   0 0.00 0.1050
#> 98   0 0.00 0.1003
#> 99   0 0.00 0.0958
#> 100  0 0.00 0.0914
#> 101  0 0.00 0.0873
#> 102  0 0.00 0.0833
#> 103  0 0.00 0.0796
#> 104  0 0.00 0.0760
#> 105  0 0.00 0.0725
#> 106  0 0.00 0.0693
#> 107  0 0.00 0.0661
#> 108  0 0.00 0.0631
#> 109  0 0.00 0.0603
#> 110  0 0.00 0.0575
#> 111  0 0.00 0.0550
#> 112  0 0.00 0.0525
#> 113  0 0.00 0.0501
#> 114  0 0.00 0.0478
#> 115  0 0.00 0.0457
#> 116  0 0.00 0.0436
#> 117  0 0.00 0.0416
#> 118  0 0.00 0.0397
#> 119  0 0.00 0.0379
#> 120  0 0.00 0.0362
#> 121  0 0.00 0.0346
#> 122  0 0.00 0.0330
#> 123  0 0.00 0.0315
#> 124  0 0.00 0.0301
#> 125  0 0.00 0.0287
#> 126  0 0.00 0.0274
#> 127  0 0.00 0.0262
#> 128  0 0.00 0.0250
#> 129  0 0.00 0.0239
#> 130  0 0.00 0.0228
#> 131  0 0.00 0.0218
#> 132  0 0.00 0.0208
#> 133  0 0.00 0.0198
#> 134  0 0.00 0.0190
#> 135  0 0.00 0.0181
#> 136  0 0.00 0.0173
#> 137  0 0.00 0.0165
#> 138  0 0.00 0.0158
#> 139  0 0.00 0.0150
#> 140  0 0.00 0.0144
#> 141  0 0.00 0.0137
#> 142  0 0.00 0.0131
#> 143  0 0.00 0.0125
#> 144  0 0.00 0.0119
#> 145  0 0.00 0.0114
#> 146  0 0.00 0.0109
#> 147  0 0.00 0.0104
#> 148  0 0.00 0.0099
#> 149  0 0.00 0.0095
#> 150  0 0.00 0.0090
#> 151  0 0.00 0.0086
#> 152  0 0.00 0.0082
#> 153  0 0.00 0.0079
#> 154  0 0.00 0.0075
#> 155  0 0.00 0.0072
#> 156  0 0.00 0.0068
#> 157  0 0.00 0.0065
#> 158  0 0.00 0.0062
#> 159  0 0.00 0.0060
#> 160  0 0.00 0.0057
#> 161  1 0.21 0.0054
#> 162  1 0.42 0.0052
#> 163  1 0.62 0.0050
#> 164  1 0.83 0.0047
#> 165  1 1.03 0.0045
#> 166  1 1.22 0.0043
#> 167  1 1.40 0.0041
#> 168  1 1.58 0.0039
#> 169  1 1.75 0.0038
#> 170  1 1.91 0.0036
#> 171  1 2.07 0.0034
#> 172  1 2.21 0.0033
#> 173  1 2.35 0.0031
#> 174  1 2.48 0.0030
#> 175  1 2.60 0.0028
#> 176  1 2.71 0.0027
#> 177  1 2.82 0.0026
#> 178  1 2.92 0.0025
#> 179  1 3.01 0.0024
#> 180  1 3.09 0.0023
#> 181  1 3.17 0.0022
#> 182  1 3.24 0.0021
#> 183  1 3.31 0.0020
#> 184  1 3.38 0.0019
#> 185  1 3.44 0.0018
#> 186  1 3.49 0.0017
#> 187  1 3.54 0.0016
#> 188  1 3.59 0.0016
#> 189  1 3.63 0.0015
#> 190  1 3.67 0.0014
#> 191  1 3.71 0.0014
#> 192  1 3.74 0.0013
#> 193  1 3.77 0.0012
#> 194  1 3.80 0.0012
#> 195  1 3.83 0.0011
#> 196  1 3.85 0.0011
#> 197  1 3.88 0.0010
#> 198  1 3.90 0.0010
#> 199  1 3.92 0.0009
#> 200  1 3.94 0.0009
#> 
#> $`0.293981481481481`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00 8785.0
#> 2    0 0.00 8388.0
#> 3    0 0.00 8009.0
#> 4    0 0.00 7646.0
#> 5    0 0.00 7300.0
#> 6    0 0.00 6970.0
#> 7    0 0.00 6655.0
#> 8    0 0.00 6354.0
#> 9    0 0.00 6067.0
#> 10   0 0.00 5792.0
#> 11   0 0.00 5530.0
#> 12   0 0.00 5280.0
#> 13   0 0.00 5041.0
#> 14   0 0.00 4813.0
#> 15   0 0.00 4596.0
#> 16   0 0.00 4388.0
#> 17   0 0.00 4189.0
#> 18   0 0.00 4000.0
#> 19   0 0.00 3819.0
#> 20   0 0.00 3646.0
#> 21   0 0.00 3481.0
#> 22   0 0.00 3324.0
#> 23   0 0.00 3174.0
#> 24   0 0.00 3030.0
#> 25   0 0.00 2893.0
#> 26   0 0.00 2762.0
#> 27   0 0.00 2637.0
#> 28   0 0.00 2518.0
#> 29   0 0.00 2404.0
#> 30   0 0.00 2295.0
#> 31   0 0.00 2191.0
#> 32   0 0.00 2092.0
#> 33   0 0.00 1998.0
#> 34   0 0.00 1907.0
#> 35   0 0.00 1821.0
#> 36   0 0.00 1739.0
#> 37   0 0.00 1660.0
#> 38   0 0.00 1585.0
#> 39   0 0.00 1513.0
#> 40   0 0.00 1445.0
#> 41   0 0.00 1380.0
#> 42   0 0.00 1317.0
#> 43   0 0.00 1258.0
#> 44   0 0.00 1201.0
#> 45   0 0.00 1146.0
#> 46   0 0.00 1095.0
#> 47   0 0.00 1045.0
#> 48   0 0.00  997.8
#> 49   0 0.00  952.6
#> 50   0 0.00  909.6
#> 51   0 0.00  868.4
#> 52   0 0.00  829.1
#> 53   0 0.00  791.6
#> 54   0 0.00  755.8
#> 55   0 0.00  721.7
#> 56   0 0.00  689.0
#> 57   0 0.00  657.8
#> 58   0 0.00  628.1
#> 59   0 0.00  599.7
#> 60   0 0.00  572.6
#> 61   0 0.00  546.7
#> 62   0 0.00  521.9
#> 63   0 0.00  498.3
#> 64   0 0.00  475.8
#> 65   0 0.00  454.3
#> 66   0 0.00  433.7
#> 67   0 0.00  414.1
#> 68   0 0.00  395.4
#> 69   0 0.00  377.5
#> 70   0 0.00  360.4
#> 71   0 0.00  344.1
#> 72   0 0.00  328.6
#> 73   0 0.00  313.7
#> 74   0 0.00  299.5
#> 75   0 0.00  286.0
#> 76   0 0.00  273.0
#> 77   0 0.00  260.7
#> 78   0 0.00  248.9
#> 79   0 0.00  237.6
#> 80   0 0.00  226.9
#> 81   0 0.00  216.6
#> 82   0 0.00  206.8
#> 83   0 0.00  197.5
#> 84   0 0.00  188.5
#> 85   0 0.00  180.0
#> 86   0 0.00  171.9
#> 87   0 0.00  164.1
#> 88   0 0.00  156.7
#> 89   0 0.00  149.6
#> 90   0 0.00  142.8
#> 91   0 0.00  136.4
#> 92   0 0.00  130.2
#> 93   0 0.00  124.3
#> 94   0 0.00  118.7
#> 95   0 0.00  113.3
#> 96   0 0.00  108.2
#> 97   0 0.00  103.3
#> 98   0 0.00   98.6
#> 99   0 0.00   94.2
#> 100  0 0.00   89.9
#> 101  0 0.00   85.8
#> 102  0 0.00   82.0
#> 103  0 0.00   78.2
#> 104  0 0.00   74.7
#> 105  0 0.00   71.3
#> 106  0 0.00   68.1
#> 107  0 0.00   65.0
#> 108  0 0.00   62.1
#> 109  0 0.00   59.3
#> 110  0 0.00   56.6
#> 111  0 0.00   54.0
#> 112  0 0.00   51.6
#> 113  0 0.00   49.3
#> 114  0 0.00   47.0
#> 115  0 0.00   44.9
#> 116  0 0.00   42.9
#> 117  0 0.00   40.9
#> 118  0 0.00   39.1
#> 119  0 0.00   37.3
#> 120  0 0.00   35.6
#> 121  0 0.00   34.0
#> 122  0 0.00   32.5
#> 123  0 0.00   31.0
#> 124  0 0.00   29.6
#> 125  0 0.00   28.3
#> 126  0 0.00   27.0
#> 127  0 0.00   25.8
#> 128  0 0.00   24.6
#> 129  0 0.00   23.5
#> 130  0 0.00   22.4
#> 131  0 0.00   21.4
#> 132  0 0.00   20.4
#> 133  0 0.00   19.5
#> 134  0 0.00   18.6
#> 135  0 0.00   17.8
#> 136  0 0.00   17.0
#> 137  0 0.00   16.2
#> 138  0 0.00   15.5
#> 139  0 0.00   14.8
#> 140  0 0.00   14.1
#> 141  0 0.00   13.5
#> 142  0 0.00   12.9
#> 143  0 0.00   12.3
#> 144  0 0.00   11.7
#> 145  0 0.00   11.2
#> 146  0 0.00   10.7
#> 147  0 0.00   10.2
#> 148  0 0.00    9.7
#> 149  0 0.00    9.3
#> 150  0 0.00    8.9
#> 151  0 0.00    8.5
#> 152  0 0.00    8.1
#> 153  0 0.00    7.7
#> 154  0 0.00    7.4
#> 155  0 0.00    7.1
#> 156  0 0.00    6.7
#> 157  0 0.00    6.4
#> 158  0 0.00    6.1
#> 159  0 0.00    5.9
#> 160  0 0.00    5.6
#> 161  0 0.00    5.3
#> 162  0 0.00    5.1
#> 163  0 0.00    4.9
#> 164  0 0.00    4.6
#> 165  0 0.00    4.4
#> 166  0 0.00    4.2
#> 167  0 0.00    4.0
#> 168  0 0.00    3.9
#> 169  0 0.00    3.7
#> 170  0 0.00    3.5
#> 171  0 0.00    3.4
#> 172  0 0.00    3.2
#> 173  0 0.00    3.1
#> 174  0 0.00    2.9
#> 175  0 0.00    2.8
#> 176  1 0.05    2.7
#> 177  1 0.26    2.5
#> 178  1 0.44    2.4
#> 179  1 0.60    2.3
#> 180  1 0.75    2.2
#> 181  1 0.88    2.1
#> 182  2 1.28    2.0
#> 183  2 1.66    1.9
#> 184  2 1.99    1.8
#> 185  2 2.28    1.8
#> 186  2 2.54    1.7
#> 187  2 2.78    1.6
#> 188  2 2.99    1.5
#> 189  2 3.17    1.5
#> 190  2 3.34    1.4
#> 191  2 3.49    1.3
#> 192  2 3.63    1.3
#> 193  2 3.75    1.2
#> 194  2 3.86    1.2
#> 195  2 3.96    1.1
#> 196  2 4.05    1.1
#> 197  2 4.13    1.0
#> 198  2 4.20    1.0
#> 199  2 4.27    0.9
#> 200  2 4.33    0.9
#> 
#> $`0.377314814814815`
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00 8776.0
#> 2    0 0.00 8379.0
#> 3    0 0.00 8000.0
#> 4    0 0.00 7639.0
#> 5    0 0.00 7293.0
#> 6    0 0.00 6963.0
#> 7    0 0.00 6648.0
#> 8    0 0.00 6348.0
#> 9    0 0.00 6061.0
#> 10   0 0.00 5786.0
#> 11   0 0.00 5525.0
#> 12   0 0.00 5275.0
#> 13   0 0.00 5036.0
#> 14   0 0.00 4808.0
#> 15   0 0.00 4591.0
#> 16   0 0.00 4383.0
#> 17   0 0.00 4185.0
#> 18   0 0.00 3996.0
#> 19   0 0.00 3815.0
#> 20   0 0.00 3643.0
#> 21   0 0.00 3478.0
#> 22   0 0.00 3321.0
#> 23   0 0.00 3170.0
#> 24   0 0.00 3027.0
#> 25   0 0.00 2890.0
#> 26   0 0.00 2759.0
#> 27   0 0.00 2635.0
#> 28   0 0.00 2515.0
#> 29   0 0.00 2402.0
#> 30   0 0.00 2293.0
#> 31   0 0.00 2189.0
#> 32   0 0.00 2090.0
#> 33   0 0.00 1996.0
#> 34   0 0.00 1905.0
#> 35   0 0.00 1819.0
#> 36   0 0.00 1737.0
#> 37   0 0.00 1658.0
#> 38   0 0.00 1583.0
#> 39   0 0.00 1512.0
#> 40   0 0.00 1443.0
#> 41   0 0.00 1378.0
#> 42   0 0.00 1316.0
#> 43   0 0.00 1256.0
#> 44   0 0.00 1199.0
#> 45   0 0.00 1145.0
#> 46   0 0.00 1093.0
#> 47   0 0.00 1044.0
#> 48   0 0.00  996.8
#> 49   0 0.00  951.7
#> 50   0 0.00  908.6
#> 51   0 0.00  867.5
#> 52   0 0.00  828.3
#> 53   0 0.00  790.8
#> 54   0 0.00  755.1
#> 55   0 0.00  720.9
#> 56   0 0.00  688.3
#> 57   0 0.00  657.2
#> 58   0 0.00  627.5
#> 59   0 0.00  599.1
#> 60   0 0.00  572.0
#> 61   0 0.00  546.1
#> 62   0 0.00  521.4
#> 63   0 0.00  497.8
#> 64   0 0.00  475.3
#> 65   0 0.00  453.8
#> 66   0 0.00  433.3
#> 67   0 0.00  413.7
#> 68   0 0.00  395.0
#> 69   0 0.00  377.1
#> 70   0 0.00  360.1
#> 71   0 0.00  343.8
#> 72   0 0.00  328.2
#> 73   0 0.00  313.4
#> 74   0 0.00  299.2
#> 75   0 0.00  285.7
#> 76   0 0.00  272.8
#> 77   0 0.00  260.4
#> 78   0 0.00  248.6
#> 79   0 0.00  237.4
#> 80   0 0.00  226.7
#> 81   0 0.00  216.4
#> 82   0 0.00  206.6
#> 83   0 0.00  197.3
#> 84   0 0.00  188.4
#> 85   0 0.00  179.8
#> 86   0 0.00  171.7
#> 87   0 0.00  163.9
#> 88   0 0.00  156.5
#> 89   0 0.00  149.4
#> 90   0 0.00  142.7
#> 91   0 0.00  136.2
#> 92   0 0.00  130.1
#> 93   0 0.00  124.2
#> 94   0 0.00  118.6
#> 95   0 0.00  113.2
#> 96   0 0.00  108.1
#> 97   0 0.00  103.2
#> 98   0 0.00   98.5
#> 99   0 0.00   94.1
#> 100  0 0.00   89.8
#> 101  0 0.00   85.8
#> 102  0 0.00   81.9
#> 103  0 0.00   78.2
#> 104  0 0.00   74.6
#> 105  0 0.00   71.3
#> 106  0 0.00   68.0
#> 107  0 0.00   65.0
#> 108  0 0.00   62.0
#> 109  0 0.00   59.2
#> 110  0 0.00   56.5
#> 111  0 0.00   54.0
#> 112  0 0.00   51.5
#> 113  0 0.00   49.2
#> 114  0 0.00   47.0
#> 115  0 0.00   44.9
#> 116  0 0.00   42.8
#> 117  0 0.00   40.9
#> 118  0 0.00   39.0
#> 119  0 0.00   37.3
#> 120  0 0.00   35.6
#> 121  0 0.00   34.0
#> 122  0 0.00   32.5
#> 123  0 0.00   31.0
#> 124  0 0.00   29.6
#> 125  0 0.00   28.2
#> 126  0 0.00   27.0
#> 127  0 0.00   25.7
#> 128  0 0.00   24.6
#> 129  0 0.00   23.5
#> 130  0 0.00   22.4
#> 131  0 0.00   21.4
#> 132  0 0.00   20.4
#> 133  0 0.00   19.5
#> 134  0 0.00   18.6
#> 135  0 0.00   17.8
#> 136  0 0.00   17.0
#> 137  0 0.00   16.2
#> 138  0 0.00   15.5
#> 139  0 0.00   14.8
#> 140  0 0.00   14.1
#> 141  0 0.00   13.5
#> 142  0 0.00   12.9
#> 143  0 0.00   12.3
#> 144  0 0.00   11.7
#> 145  0 0.00   11.2
#> 146  0 0.00   10.7
#> 147  0 0.00   10.2
#> 148  0 0.00    9.7
#> 149  0 0.00    9.3
#> 150  0 0.00    8.9
#> 151  0 0.00    8.5
#> 152  0 0.00    8.1
#> 153  0 0.00    7.7
#> 154  0 0.00    7.4
#> 155  0 0.00    7.0
#> 156  0 0.00    6.7
#> 157  0 0.00    6.4
#> 158  0 0.00    6.1
#> 159  0 0.00    5.9
#> 160  0 0.00    5.6
#> 161  1 0.21    5.3
#> 162  1 0.42    5.1
#> 163  1 0.62    4.9
#> 164  1 0.83    4.6
#> 165  1 1.03    4.4
#> 166  1 1.22    4.2
#> 167  1 1.40    4.0
#> 168  1 1.58    3.9
#> 169  1 1.75    3.7
#> 170  1 1.91    3.5
#> 171  1 2.07    3.4
#> 172  1 2.21    3.2
#> 173  1 2.35    3.1
#> 174  1 2.48    2.9
#> 175  1 2.60    2.8
#> 176  1 2.71    2.7
#> 177  1 2.82    2.5
#> 178  2 2.95    2.4
#> 179  2 3.38    2.3
#> 180  2 3.76    2.2
#> 181  2 4.11    2.1
#> 182  2 4.42    2.0
#> 183  2 4.71    1.9
#> 184  2 4.96    1.8
#> 185  2 5.19    1.8
#> 186  2 5.40    1.7
#> 187  2 5.59    1.6
#> 188  2 5.76    1.5
#> 189  2 5.92    1.5
#> 190  2 6.06    1.4
#> 191  2 6.18    1.3
#> 192  2 6.30    1.3
#> 193  2 6.41    1.2
#> 194  2 6.50    1.2
#> 195  2 6.59    1.1
#> 196  2 6.67    1.1
#> 197  2 6.74    1.0
#> 198  2 6.81    1.0
#> 199  2 6.87    0.9
#> 200  2 6.92    0.9
#> 
esm_max_t1$predictors
#> # A tibble: 1 × 8
#>   c1    c2    c3    c4      c5      c6    c7    c8   
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr>
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth
esm_max_t1$performance
#> # A tibble: 7 × 33
#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr>   <chr>            <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 esm_max equal_sens_…     0.555          10         10        1      0        1
#> 2 esm_max lpt              0.551          10         10        1      0        1
#> 3 esm_max max_fpb          0.551          10         10        1      0        1
#> 4 esm_max max_jaccard      0.551          10         10        1      0        1
#> 5 esm_max max_sens_sp…     0.581          10         10        1      0        1
#> 6 esm_max max_sorensen     0.551          10         10        1      0        1
#> 7 esm_max sensitivity      0.592          10         10        1      0        1
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
esm_max_t1$performance_part
#> # A tibble: 42 × 21
#>    model replicates part  threshold thr_value n_presences n_absences   TPR   TNR
#>    <chr> <chr>      <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#>  1 esm_… .part1     1     max_sore…     0.551           4          4     1     1
#>  2 esm_… .part1     1     max_jacc…     0.551           4          4     1     1
#>  3 esm_… .part1     1     max_fpb       0.551           4          4     1     1
#>  4 esm_… .part1     1     max_sens…     0.551           4          4     1     1
#>  5 esm_… .part1     1     equal_se…     0.551           4          4     1     1
#>  6 esm_… .part1     1     lpt           0.551           4          4     1     1
#>  7 esm_… .part1     1     sensitiv…     0.551           4          4     1     1
#>  8 esm_… .part1     2     max_sore…     0.715           3          3     1     1
#>  9 esm_… .part1     2     max_jacc…     0.715           3          3     1     1
#> 10 esm_… .part1     2     max_fpb       0.715           3          3     1     1
#> # ℹ 32 more rows
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
# }
```
