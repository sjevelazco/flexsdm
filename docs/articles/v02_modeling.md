# flexsdm: Overview of Modeling functions

## Introduction

Species distribution modeling (SDM) has become a standard tool in
multiple research areas, including ecology, conservation biology,
biogeography, paleobiogeography, and epidemiology. SDM is an area of
active theoretical and methodological research. The *flexsdm* package
provides users the ability to manipulate and parameterize models in a
variety of ways that meet their unique research needs.

This flexibility enables users to define their own complete or partial
modeling procedure specific for their modeling situations (e.g., number
of variables, number of records, different algorithms and ensemble
methods, algorithms tuning, etc.).

In this vignette, users will learn about the second set of functions in
the *flexsdm* package that fall under the “modeling” umbrella. These
functions were designed to construct and validate different types of
models and can be grouped into fit\_\* , tune\_\* , and esm\_\* family
functions. In addition there is a function to perform ensemble modeling.

The fit\_\* functions construct and validate models with default
hyper-parameter values. The tune\_\* functions construct and validate
models by searching for the best combination of hyper-parameter values,
and esm\_ functions can be used for constructing and validating Ensemble
of Small Models. Finally, the fit_ensemble() function is for fitting and
validating ensemble models.

These are the functions for model construction and validation:

**fit\_\* functions family**

- fit_gam() Fit and validate Generalized Additive Models

- fit_gau() Fit and validate Gaussian Process models

- fit_gbm() Fit and validate Generalized Boosted Regression models

- fit_glm() Fit and validate Generalized Linear Models

- fit_max() Fit and validate Maximum Entropy models

- fit_net() Fit and validate Neural Networks models

- fit_raf() Fit and validate Random Forest models

- fit_svm() Fit and validate Support Vector Machine models

**tune\_\* functions family**

- tune_gbm() Fit and validate Generalized Boosted Regression models with
  exploration of hyper-parameters

- tune_max() Fit and validate Maximum Entropy models with exploration of
  hyper-parameters

- tune_net() Fit and validate Neural Networks models with exploration of
  hyper-parameters

- tune_raf() Fit and validate Random Forest models with exploration of
  hyper-parameters

- tune_svm() Fit and validate Support Vector Machine models with
  exploration of hyper-parameters

**model ensemble**

- fit_ensemble() Fit and validate ensemble models with different
  ensemble methods

**esm\_\* functions family**

- esm_gam() Fit and validate Generalized Additive Models with Ensemble
  of Small Model approach

- esm_gau() Fit and validate Gaussian Process models Models with
  Ensemble of Small Model approach

- esm_gbm() Fit and validate Generalized Boosted Regression models with
  Ensemble of Small Model approach

- esm_glm() Fit and validate Generalized Linear Models with Ensemble of
  Small Model approach

- esm_max() Fit and validate Maximum Entropy models with Ensemble of
  Small Model approach

- esm_net() Fit and validate Neural Networks models with Ensemble of
  Small Model approach

- esm_svm() Fit and validate Support Vector Machine models with Ensemble
  of Small Model approach

## Installation

First, install the flexsdm package. You can install the released version
of *flexsdm* from [github](https://github.com/sjevelazco/flexsdm) with:

``` r
# devtools::install_github('sjevelazco/flexsdm')
require(flexsdm)
#> Loading required package: flexsdm
require(terra)
#> Loading required package: terra
#> terra 1.8.60
#> 
#> Attaching package: 'terra'
#> The following object is masked from 'package:knitr':
#> 
#>     spin
require(dplyr)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:terra':
#> 
#>     intersect, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Project directory setup

Decide where on your computer you would like to store the inputs and
outputs of your project (this will be your main directory). Use an
existing one or use dir.create() to create your main directory. Then
specify whether or not to include folders for projections, calibration
areas, algorithms, ensembles, and thresholds. For more details see
[Vignette
01_pre_modeling](https://sjevelazco.github.io/flexsdm/articles/01_pre_modeling.html)

## Data, species occurrence and background data

In this tutorial, we will be using species occurrences and environmental
data that are available through the *flexsdm* package. The “abies”
example dataset includes a pr_ab column (presence = 1, and absence = 0),
location columns (x, y) and other environmental data. You can load the
“abies” data into your local R environment by using the code below:

(THIS EXAMPLE LOOKS A LITTLE STRANGE BECAUSE WE ARE ALSO USING
BACKGROUND DATA, WHILE THE ABIES DATASET CLEARLY HAS ABSENCES…)

``` r
data("abies")
data("backg")

dplyr::glimpse(abies)
#> Rows: 1,400
#> Columns: 13
#> $ id       <int> 715, 5680, 7907, 1850, 1702, 10036, 12384, 6513, 9884, 8651, …
#> $ pr_ab    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ x        <dbl> -95417.134, 98986.536, 121474.257, -39976.221, 111372.261, -2…
#> $ y        <dbl> 314240.13, -159415.18, -99463.44, -17456.11, -91404.05, 39222…
#> $ aet      <dbl> 323.1133, 447.5567, 182.2833, 372.3867, 209.4567, 308.3000, 5…
#> $ cwd      <dbl> 546.1400, 815.4033, 271.1800, 946.2933, 398.5500, 534.9533, 3…
#> $ tmin     <dbl> 1.2433, 9.4267, -4.9500, 8.7767, -4.0333, 4.6600, 4.3800, 4.9…
#> $ ppt_djf  <dbl> 62.7257, 129.6406, 150.7003, 116.0236, 164.9327, 166.2220, 48…
#> $ ppt_jja  <dbl> 17.7941, 6.4317, 11.2294, 2.7020, 9.2686, 16.5310, 41.2494, 8…
#> $ pH       <dbl> 5.773341, 5.600000, 0.000000, 6.411796, 0.000000, 5.700000, 5…
#> $ awc      <dbl> 0.10837019, 0.16000000, 0.00000000, 0.09719457, 0.00000000, 0…
#> $ depth    <dbl> 152.000000, 201.000000, 0.000000, 59.759930, 0.000000, 112.99…
#> $ landform <fct> 7, 11, 15, 14, 15, 15, 7, 15, 4, 10, 6, 10, 10, 15, 10, 11, 1…
dplyr::glimpse(backg)
#> Rows: 5,000
#> Columns: 13
#> $ pr_ab        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
#> $ x            <dbl> 160779.16, 36849.16, -240170.84, -152420.84, -193190.84, …
#> $ y            <dbl> -449968.33, 24151.67, 90031.67, -143518.33, 24151.67, 223…
#> $ aet          <dbl> 280.4567, 259.7800, 400.1767, 367.4833, 397.3667, 385.263…
#> $ cwd          <dbl> 1137.2433, 381.5367, 699.6500, 843.4467, 842.3833, 637.35…
#> $ tmin         <dbl> 13.5100, -3.1733, 8.6800, 9.0133, 8.9700, 4.9333, 6.2933,…
#> $ ppt_djf      <dbl> 71.2741, 171.4537, 285.0893, 72.0309, 125.2467, 226.1534,…
#> $ ppt_jja      <dbl> 1.1920, 17.5193, 5.0158, 1.2047, 1.9778, 8.1554, 18.4182,…
#> $ pH           <dbl> 0.0000000, 0.2122687, 5.7222223, 7.5350823, 6.1963525, 5.…
#> $ awc          <dbl> 0.000000000, 0.003473487, 0.080370426, 0.170000002, 0.131…
#> $ depth        <dbl> 0.00000, 201.00000, 50.07409, 154.39426, 122.39575, 56.17…
#> $ percent_clay <dbl> 0.0000000, 0.4438345, 18.4111176, 46.9751244, 37.1873169,…
#> $ landform     <fct> 13, 10, 6, 6, 10, 14, 8, 14, 6, 7, 11, 14, 14, 10, 6, 6, …
```

If you want to replace the abies dataset with your own data, make sure
that your dataset contains the environmental conditions related to
presence-absence data. We use a pre-modeling family function for a
k-fold partition method (to be used for cross-validation). In the
partition method the number of folds or replications must be the same
for presence-absence and for background points datasets.

``` r
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

dplyr::glimpse(abies2)
#> Rows: 1,400
#> Columns: 14
#> $ id       <int> 715, 5680, 7907, 1850, 1702, 10036, 12384, 6513, 9884, 8651, …
#> $ pr_ab    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ x        <dbl> -95417.134, 98986.536, 121474.257, -39976.221, 111372.261, -2…
#> $ y        <dbl> 314240.13, -159415.18, -99463.44, -17456.11, -91404.05, 39222…
#> $ aet      <dbl> 323.1133, 447.5567, 182.2833, 372.3867, 209.4567, 308.3000, 5…
#> $ cwd      <dbl> 546.1400, 815.4033, 271.1800, 946.2933, 398.5500, 534.9533, 3…
#> $ tmin     <dbl> 1.2433, 9.4267, -4.9500, 8.7767, -4.0333, 4.6600, 4.3800, 4.9…
#> $ ppt_djf  <dbl> 62.7257, 129.6406, 150.7003, 116.0236, 164.9327, 166.2220, 48…
#> $ ppt_jja  <dbl> 17.7941, 6.4317, 11.2294, 2.7020, 9.2686, 16.5310, 41.2494, 8…
#> $ pH       <dbl> 5.773341, 5.600000, 0.000000, 6.411796, 0.000000, 5.700000, 5…
#> $ awc      <dbl> 0.10837019, 0.16000000, 0.00000000, 0.09719457, 0.00000000, 0…
#> $ depth    <dbl> 152.000000, 201.000000, 0.000000, 59.759930, 0.000000, 112.99…
#> $ landform <fct> 7, 11, 15, 14, 15, 15, 7, 15, 4, 10, 6, 10, 10, 15, 10, 11, 1…
#> $ .part    <int> 3, 2, 5, 3, 5, 4, 2, 1, 2, 2, 1, 2, 1, 2, 2, 2, 5, 4, 5, 1, 5…
```

Now, in the abies2 object we have a new column called “.part” with the 5
k-folds (1, 2, 3, 4, 5), indicating which partition each record (row) is
a member of. Next, we have to apply the same partition method and number
of folds to the environmental conditions of the background points.

``` r
backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

dplyr::glimpse(backg2)
#> Rows: 5,000
#> Columns: 14
#> $ pr_ab        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
#> $ x            <dbl> 160779.16, 36849.16, -240170.84, -152420.84, -193190.84, …
#> $ y            <dbl> -449968.33, 24151.67, 90031.67, -143518.33, 24151.67, 223…
#> $ aet          <dbl> 280.4567, 259.7800, 400.1767, 367.4833, 397.3667, 385.263…
#> $ cwd          <dbl> 1137.2433, 381.5367, 699.6500, 843.4467, 842.3833, 637.35…
#> $ tmin         <dbl> 13.5100, -3.1733, 8.6800, 9.0133, 8.9700, 4.9333, 6.2933,…
#> $ ppt_djf      <dbl> 71.2741, 171.4537, 285.0893, 72.0309, 125.2467, 226.1534,…
#> $ ppt_jja      <dbl> 1.1920, 17.5193, 5.0158, 1.2047, 1.9778, 8.1554, 18.4182,…
#> $ pH           <dbl> 0.0000000, 0.2122687, 5.7222223, 7.5350823, 6.1963525, 5.…
#> $ awc          <dbl> 0.000000000, 0.003473487, 0.080370426, 0.170000002, 0.131…
#> $ depth        <dbl> 0.00000, 201.00000, 50.07409, 154.39426, 122.39575, 56.17…
#> $ percent_clay <dbl> 0.0000000, 0.4438345, 18.4111176, 46.9751244, 37.1873169,…
#> $ landform     <fct> 13, 10, 6, 6, 10, 14, 8, 14, 6, 7, 11, 14, 14, 10, 6, 6, …
#> $ .part        <int> 4, 4, 1, 5, 5, 2, 5, 3, 2, 5, 4, 1, 4, 1, 5, 1, 1, 5, 4, …
```

In backg2 object we have a new column called “.part” with the 5 k-folds
(1, 2, 3, 4, 5).

### 1. Fit and validate models

We fit and validate models: I. a maximum entropy model with default
hyper-parameter values (flexsdm::fit_max) and II. a random forest model
with exploration of hyper-parameters (flexsdm::tune_raf).

I. Maximum Entropy models with default hyper-parameter values.

``` r
max_t1 <- fit_max(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  background = backg2,
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
  clamp = TRUE,
  classes = "default",
  pred_type = "cloglog",
  regmult = 1
)
#> Formula used for model fitting:
#> ~aet + ppt_jja + pH + awc + depth + I(aet^2) + I(ppt_jja^2) + I(pH^2) + I(awc^2) + I(depth^2) + hinge(aet) + hinge(ppt_jja) + hinge(pH) + hinge(awc) + hinge(depth) + ppt_jja:aet + pH:aet + awc:aet + depth:aet + pH:ppt_jja + awc:ppt_jja + depth:ppt_jja + awc:pH + depth:pH + depth:awc + categorical(landform) - 1
#> Replica number: 1/1
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
```

This function returns a list object with the following elements:

``` r
names(max_t1)
#> [1] "model"            "predictors"       "performance"      "performance_part"
#> [5] "data_ens"
```

model: A “MaxEnt” class object. This object can be used for predicting.

``` r
options(max.print = 20)
max_t1$model
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df  %Dev  Lambda
#> 1    0  0.00 21.3700
#> 2    0  0.00 20.4100
#> 3    0  0.00 19.4800
#> 4    0  0.00 18.6000
#> 5    0  0.00 17.7600
#> 6    0  0.00 16.9600
#>  [ reached getOption("max.print") -- omitted 194 rows ]
```

predictors: A tibble with quantitative (c column names) and qualitative
(f column names) variables use for modeling.

``` r
max_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
```

performance: The performance metric (see sdm_eval). Those metrics that
are threshold dependent are calculated based on the threshold specified
in the argument. We can see all the selected threshold values.

``` r
max_t1$performance
#> # A tibble: 3 × 25
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 max   equal_sens_sp…     0.573         700        700    0.674 0.0164    0.674
#> 2 max   max_sens_spec      0.416         700        700    0.909 0.0260    0.52 
#> 3 max   max_sorensen       0.335         700        700    0.95  0.0101    0.469
#> # ℹ 17 more variables: TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,
#> #   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,
#> #   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, AUC_mean <dbl>,
#> #   AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>, IMAE_mean <dbl>,
#> #   IMAE_sd <dbl>
```

Predicted suitability for each test partition (row) based on the best
model. This database is used in fit_ensemble.

``` r
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab   pred
#>    <chr>  <chr>      <chr> <dbl>  <dbl>
#>  1 8      .part      1         0 0.600 
#>  2 11     .part      1         0 0.237 
#>  3 13     .part      1         0 0.0483
#>  4 20     .part      1         0 0.115 
#>  5 32     .part      1         0 0.716 
#>  6 33     .part      1         0 0.0430
#>  7 48     .part      1         0 0.143 
#>  8 55     .part      1         0 0.726 
#>  9 65     .part      1         0 0.850 
#> 10 75     .part      1         0 0.308 
#> # ℹ 1,390 more rows
```

II- Random forest models with exploration of hyper-parameters.

First, we create a data.frame that provides hyper-parameters values to
be tested. It Is recommended to generate this data.frame.
Hyper-parameter needed for tuning is ‘mtry’. The maximum mtry must be
equal to total number of predictors.

``` r
tune_grid <-
  expand.grid(
    mtry = seq(1, 7, 1),
    ntree = c(300, 500, 700, 900)
  )
```

We use the same data object abies2, with the same k-fold partition
method:

``` r
rf_t <-
  tune_raf(
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
```

Let’s see what the output object contains. This function returns a list
object with the following elements:

``` r
names(rf_t)
#> [1] "model"             "predictors"        "performance"      
#> [4] "performance_part"  "hyper_performance" "data_ens"
```

model: A “randomForest” class object. This object can be used to see the
formula details, a basic summary o fthe model, and for predicting.

``` r
rf_t$model
#> 
#> Call:
#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = 500,      importance = TRUE, ) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 2
#> 
#>         OOB estimate of  error rate: 10.93%
#> Confusion matrix:
#>     0   1 class.error
#> 0 606  94  0.13428571
#> 1  59 641  0.08428571
```

predictors: A tibble with quantitative (c column names) and qualitative
(f column names) variables use for modeling.

``` r
rf_t$predictors
#> # A tibble: 1 × 9
#>   c1    c2    c3    c4      c5      c6    c7    c8    f       
#>   <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth landform
```

performance: The performance metric (see sdm_eval). Those metrics that
are threshold dependent are calculated based on the threshold specified
in the argument. We can see all the selected threshold values.

``` r
rf_t$performance
#> # A tibble: 1 × 27
#>    mtry ntree model threshold   thr_value n_presences n_absences TPR_mean TPR_sd
#>   <dbl> <dbl> <chr> <chr>           <dbl>       <int>      <int>    <dbl>  <dbl>
#> 1     2   300 raf   max_sens_s…      0.53         700        700    0.913 0.0383
#> # ℹ 18 more variables: TNR_mean <dbl>, TNR_sd <dbl>, SORENSEN_mean <dbl>,
#> #   SORENSEN_sd <dbl>, JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>,
#> #   FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>,
#> #   AUC_mean <dbl>, AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
```

Predicted suitability for each test partition (row) based on the best
model. This database is used in fit_ensemble.

``` r
rf_t$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <fct> <dbl>
#>  1 8      .part      1     0     0.156
#>  2 11     .part      1     0     0.152
#>  3 13     .part      1     0     0.01 
#>  4 20     .part      1     0     0.476
#>  5 32     .part      1     0     0.126
#>  6 33     .part      1     0     0.078
#>  7 48     .part      1     0     0.01 
#>  8 55     .part      1     0     0.148
#>  9 65     .part      1     0     0.444
#> 10 75     .part      1     0     0.086
#> # ℹ 1,390 more rows
```

These model objects can be used in flexsdm::fit_ensemble().

### 2. Model Ensemble

In this example we fit and validate and ensemble model using the two
model objects that were just created.

``` r
# Fit and validate ensemble model
an_ensemble <- fit_ensemble(
  models = list(max_t1, rf_t),
  ens_method = "meansup",
  thr = NULL,
  thr_model = "max_sens_spec",
  metric = "TSS"
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
```

``` r
# Outputs
names(an_ensemble)
#> [1] "models"           "thr_metric"       "predictors"       "performance"     
#> [5] "performance_part"

an_ensemble$thr_metric
#> [1] "max_sens_spec" "TSS_mean"
an_ensemble$predictors
#> # A tibble: 2 × 9
#>   c1    c2      c3    c4      c5      f        c6    c7    c8   
#>   <chr> <chr>   <chr> <chr>   <chr>   <chr>    <chr> <chr> <chr>
#> 1 aet   ppt_jja pH    awc     depth   landform NA    NA    NA   
#> 2 aet   cwd     tmin  ppt_djf ppt_jja landform pH    awc   depth
an_ensemble$performance
#> # A tibble: 7 × 25
#>   model   threshold   thr_value n_presences n_absences TPR_mean  TPR_sd TNR_mean
#>   <chr>   <chr>           <dbl>       <int>      <int>    <dbl>   <dbl>    <dbl>
#> 1 meansup equal_sens…     0.58          700        700    0.876 0.0251     0.876
#> 2 meansup lpt             0.028         700        700    1     0          0.421
#> 3 meansup max_fpb         0.48          700        700    0.933 0.0310     0.843
#> 4 meansup max_jaccard     0.48          700        700    0.933 0.0310     0.843
#> 5 meansup max_sens_s…     0.48          700        700    0.91  0.0345     0.87 
#> 6 meansup max_sorens…     0.48          700        700    0.933 0.0310     0.843
#> 7 meansup sensitivity     0.534         700        700    0.897 0.00391    0.857
#> # ℹ 17 more variables: TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,
#> #   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,
#> #   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, AUC_mean <dbl>,
#> #   AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>, IMAE_mean <dbl>,
#> #   IMAE_sd <dbl>
```

### 3. Fit and validate models with Ensemble of Small Model approach

This method consists of creating bivariate models with all pair-wise
combinations of predictors and perform an ensemble based on the average
of suitability weighted by Somers’ D metric (D = 2 x (AUC -0.5)). ESM is
recommended for modeling species with very few occurrences. This
function does not allow categorical variables because the use of these
types of variables could be problematic when applied to species with few
occurrences. For more detail see Breiner et al. (2015, 2018)

``` r
data("abies")
library(dplyr)

# Create a smaller subset of occurrences
set.seed(10)
abies2 <- abies %>%
  na.omit() %>%
  group_by(pr_ab) %>%
  dplyr::slice_sample(n = 10) %>%
  group_by()
```

We can use different methods in the flexsdm::part_random function
according to our data. See
[part_random](https://sjevelazco.github.io/flexsdm/reference/part_random.html)
for more details.

``` r
# Using k-fold partition method for model cross validation
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
```

This function constructs Generalized Additive Models using the Ensembles
of Small Models (ESM) approach (Breiner et al., 2015, 2018).

``` r
# We set the model without threshold specification and with the kfold created above
esm_gam_t1 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#> 
#> Model has more coefficients than data used for training it. Try to reduce k
```

This function returns a list object with the following elements:

``` r
names(esm_gam_t1)
#> NULL
```

esm_model: A list with “GAM” class object for each bivariate model. This
object can be used for predicting using the ESM approachwith sdm_predict
function.

``` r
options(max.print = 10) # If you don't want to see printed all the output
esm_gam_t1$esm_model
#> NULL
```

predictors: A tibble with variables use for modeling.

``` r
esm_gam_t1$predictors
#> NULL
```

performance: Performance metric (see sdm_eval). Those threshold
dependent metrics are calculated based on the threshold specified in the
argument.

``` r
esm_gam_t1$performance
#> NULL
```

Now, we test the rep_kfold partition method. In this method ‘folds’
refers to the number of partitions for data partitioning and ‘replicate’
refers to the number of replicates. Both assume values \>=1.

``` r
# Remove the previous k-fold partition
abies2 <- abies2 %>% select(-starts_with("."))

# Test with rep_kfold partition using 3 folds and 5 replicates
set.seed(10)
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
```

We use the new rep_kfold partition in the gam model

``` r
esm_gam_t2 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#> 
#> Model has more coefficients than data used for training it. Try to reduce k
```

Test with random bootstrap partitioning. In method ‘replicate’ refers to
the number of replicates (assumes a value \>=1), ‘proportion’ refers to
the proportion of occurrences used for model fitting (assumes a value
\>0 and \<=1). With this method we can configure the proportion of
training and testing data according to the species occurrences. In this
example, proportion=‘0.7’ indicates that 70% of data will be used for
model training, while 30% will be used for model testing. For this
method, the function will return .partX columns with “train” or “test”
words as the entries.

``` r
# Remove the previous k-fold partition
abies2 <- abies2 %>% select(-starts_with("."))

# Test with bootstrap partition using 10 replicates
set.seed(10)
abies2 <- part_random(
  data = abies2,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 10, proportion = 0.7)
)
abies2
#> # A tibble: 20 × 23
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
#> # ℹ 12 more variables: depth <dbl>, landform <fct>, .part1 <chr>, .part2 <chr>,
#> #   .part3 <chr>, .part4 <chr>, .part5 <chr>, .part6 <chr>, .part7 <chr>,
#> #   .part8 <chr>, .part9 <chr>, .part10 <chr>
```

Use the new rep_kfold partition in the gam model

``` r
esm_gam_t3 <- esm_gam(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
  partition = ".part",
  thr = NULL
)
#> 
#> Model has more coefficients than data used for training it. Try to reduce k
```

\#=========#=========#=========#=========#=========#=========#=========#

**Vignette still under construction and changes**

\#=========#=========#=========#=========#=========#=========#=========#
