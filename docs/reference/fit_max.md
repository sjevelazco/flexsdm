# Fit and validate Maximum Entropy models

Fit and validate Maximum Entropy models

## Usage

``` r
fit_max(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition = NULL,
  background = NULL,
  thr = NULL,
  clamp = TRUE,
  classes = "default",
  pred_type = "cloglog",
  regmult = 1
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

  formula. A formula object with response and predictor variables. See
  maxnet.formula function from maxnet package. Note that the variables
  used here must be consistent with those used in response, predictors,
  and predictors_f arguments. Default NULL.

- partition:

  character. Column name with training and validation partition groups.
  If partition = NULL, the model will be validated with the same data
  used for fitting.

- background:

  data.frame. Database including only those rows with 0 values in the
  response column and the predictors variables. All column names must be
  consistent with data. Default NULL

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

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If a sensitivity values is not
    specified the default used is 0.9.

  If more than one threshold type is used they must be concatenated,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- clamp:

  logical. If TRUE, predictors and features are restricted to the range
  seen during model training.

- classes:

  character. A single feature of any combinations of them. Features are
  symbolized by letters: l (linear), q (quadratic), h (hinge), p
  (product), and t (threshold). Usage classes = "lpq". Default "default"
  (see details).

- pred_type:

  character. Type of response required available "link", "exponential",
  "cloglog" and "logistic". Default "cloglog"

- regmult:

  numeric. A constant to adjust regularization. Default 1.

## Value

A list object with:

- model: A "maxnet" class object from maxnet package. This object can be
  used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).
  Threshold dependent metrics are calculated based on the threshold
  specified in thr argument.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- data_ens: Predicted suitability for each test partition based on the
  best model. This database is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## Details

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

## See also

[`fit_gam`](https://sjevelazco.github.io/flexsdm/reference/fit_gam.md),
[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md),
[`fit_gbm`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md),
[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md),
[`fit_net`](https://sjevelazco.github.io/flexsdm/reference/fit_net.md),
[`fit_raf`](https://sjevelazco.github.io/flexsdm/reference/fit_raf.md),
and
[`fit_svm`](https://sjevelazco.github.io/flexsdm/reference/fit_svm.md).

## Examples

``` r
# \donttest{
data("abies")
data("backg")
set.seed(1)
backg <- backg[sample(nrow(backg), 1000), ] # subsample to speed up this example
abies # environmental conditions of presence-absence data
#> # A tibble: 1,400 × 13
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
#> # ℹ 2 more variables: depth <dbl>, landform <fct>
backg # environmental conditions of background points
#> # A tibble: 1,000 × 13
#>    pr_ab        x        y   aet   cwd  tmin ppt_djf ppt_jja    pH     awc depth
#>    <dbl>    <dbl>    <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl>   <dbl> <dbl>
#>  1     0   23889. -320098.  194. 1184.  7.00    36.7    1.01 8.12  0.167   201  
#>  2     0  278769. -439708.  369. 1049.  7.70   101.     8.18 6.60  0.120    36  
#>  3     0  110019. -208858.  370. 1014. 10.0     88.4    2.90 5.65  0.0727  123. 
#>  4     0 -163491.  213962.  342.  896.  9.59   127.     5.94 6.40  0.0900   46  
#>  5     0 -217491.  230702.  331.  907. 10.8    130.     5.91 6     0.120   201  
#>  6     0 -262581.  287402.  306.  621.  2.86   169.     8.63 6.32  0.0936  102. 
#>  7     0 -191841.  289022.  397.  780.  9.95   184.    10.2  5.30  0.0913   67.3
#>  8     0  107049. -324958.  211. 1314. 12.2     34.0    1.20 7.80  0.140   201  
#>  9     0   -7701.   70592.  284.  397. -1.21   238.    14.5  0.992 0.00740 173. 
#> 10     0 -245841.  391622.  246.  663.  3.10   116.     9.60 5.89  0.0800  118. 
#> # ℹ 990 more rows
#> # ℹ 2 more variables: percent_clay <dbl>, landform <fct>

# Using k-fold partition method
# Note that the partition method, number of folds or replications must
# be the same for presence-absence and background points datasets
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
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

backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
backg2
#> # A tibble: 1,000 × 14
#>    pr_ab        x        y   aet   cwd  tmin ppt_djf ppt_jja    pH     awc depth
#>    <dbl>    <dbl>    <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl>   <dbl> <dbl>
#>  1     0   23889. -320098.  194. 1184.  7.00    36.7    1.01 8.12  0.167   201  
#>  2     0  278769. -439708.  369. 1049.  7.70   101.     8.18 6.60  0.120    36  
#>  3     0  110019. -208858.  370. 1014. 10.0     88.4    2.90 5.65  0.0727  123. 
#>  4     0 -163491.  213962.  342.  896.  9.59   127.     5.94 6.40  0.0900   46  
#>  5     0 -217491.  230702.  331.  907. 10.8    130.     5.91 6     0.120   201  
#>  6     0 -262581.  287402.  306.  621.  2.86   169.     8.63 6.32  0.0936  102. 
#>  7     0 -191841.  289022.  397.  780.  9.95   184.    10.2  5.30  0.0913   67.3
#>  8     0  107049. -324958.  211. 1314. 12.2     34.0    1.20 7.80  0.140   201  
#>  9     0   -7701.   70592.  284.  397. -1.21   238.    14.5  0.992 0.00740 173. 
#> 10     0 -245841.  391622.  246.  663.  3.10   116.     9.60 5.89  0.0800  118. 
#> # ℹ 990 more rows
#> # ℹ 3 more variables: percent_clay <dbl>, landform <fct>, .part <int>

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
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
length(max_t1)
#> [1] 5

max_t1$model
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev Lambda
#> 1    0 0.00 69.040
#> 2    0 0.00 65.910
#> 3    0 0.00 62.930
#> 4    0 0.00 60.090
#> 5    0 0.00 57.370
#> 6    0 0.00 54.770
#> 7    0 0.00 52.300
#> 8    0 0.00 49.930
#> 9    0 0.00 47.670
#> 10   0 0.00 45.520
#> 11   0 0.00 43.460
#> 12   0 0.00 41.490
#> 13   0 0.00 39.620
#> 14   0 0.00 37.820
#> 15   0 0.00 36.110
#> 16   0 0.00 34.480
#> 17   0 0.00 32.920
#> 18   0 0.00 31.430
#> 19   0 0.00 30.010
#> 20   0 0.00 28.650
#> 21   0 0.00 27.360
#> 22   0 0.00 26.120
#> 23   0 0.00 24.940
#> 24   0 0.00 23.810
#> 25   0 0.00 22.730
#> 26   0 0.00 21.710
#> 27   0 0.00 20.720
#> 28   0 0.00 19.790
#> 29   0 0.00 18.890
#> 30   0 0.00 18.040
#> 31   0 0.00 17.220
#> 32   0 0.00 16.440
#> 33   0 0.00 15.700
#> 34   0 0.00 14.990
#> 35   0 0.00 14.310
#> 36   0 0.00 13.660
#> 37   0 0.00 13.050
#> 38   0 0.00 12.460
#> 39   0 0.00 11.890
#> 40   0 0.00 11.350
#> 41   0 0.00 10.840
#> 42   0 0.00 10.350
#> 43   0 0.00  9.882
#> 44   0 0.00  9.435
#> 45   0 0.00  9.008
#> 46   0 0.00  8.601
#> 47   0 0.00  8.212
#> 48   0 0.00  7.841
#> 49   0 0.00  7.486
#> 50   0 0.00  7.147
#> 51   0 0.00  6.824
#> 52   0 0.00  6.516
#> 53   0 0.00  6.221
#> 54   0 0.00  5.939
#> 55   0 0.00  5.671
#> 56   0 0.00  5.414
#> 57   0 0.00  5.169
#> 58   0 0.00  4.936
#> 59   0 0.00  4.712
#> 60   0 0.00  4.499
#> 61   0 0.00  4.296
#> 62   0 0.00  4.102
#> 63   0 0.00  3.916
#> 64   0 0.00  3.739
#> 65   0 0.00  3.570
#> 66   0 0.00  3.408
#> 67   0 0.00  3.254
#> 68   0 0.00  3.107
#> 69   0 0.00  2.966
#> 70   0 0.00  2.832
#> 71   0 0.00  2.704
#> 72   0 0.00  2.582
#> 73   0 0.00  2.465
#> 74   0 0.00  2.354
#> 75   0 0.00  2.247
#> 76   0 0.00  2.146
#> 77   0 0.00  2.049
#> 78   0 0.00  1.956
#> 79   0 0.00  1.867
#> 80   0 0.00  1.783
#> 81   0 0.00  1.702
#> 82   0 0.00  1.625
#> 83   0 0.00  1.552
#> 84   1 0.21  1.482
#> 85   1 0.48  1.415
#> 86   1 0.73  1.351
#> 87   1 0.97  1.290
#> 88   1 1.19  1.231
#> 89   1 1.41  1.176
#> 90   1 1.61  1.122
#> 91   1 1.80  1.072
#> 92   1 1.98  1.023
#> 93   1 2.15  0.977
#> 94   1 2.31  0.933
#> 95   1 2.46  0.890
#> 96   1 2.61  0.850
#> 97   1 2.75  0.812
#> 98   1 2.88  0.775
#> 99   1 3.00  0.740
#> 100  1 3.12  0.707
#> 101  1 3.23  0.675
#> 102  1 3.34  0.644
#> 103  1 3.44  0.615
#> 104  1 3.54  0.587
#> 105  1 3.63  0.561
#> 106  1 3.72  0.535
#> 107  1 3.80  0.511
#> 108  1 3.88  0.488
#> 109  1 3.95  0.466
#> 110  1 4.03  0.445
#> 111  2 4.12  0.425
#> 112  2 4.21  0.405
#> 113  2 4.29  0.387
#> 114  2 4.37  0.370
#> 115  2 4.45  0.353
#> 116  2 4.52  0.337
#> 117  2 4.58  0.322
#> 118  3 4.67  0.307
#> 119  3 4.75  0.293
#> 120  3 4.82  0.280
#> 121  4 4.90  0.267
#> 122  4 5.01  0.255
#> 123  3 5.10  0.244
#> 124  3 5.18  0.233
#> 125  3 5.26  0.222
#> 126  3 5.33  0.212
#> 127  3 5.40  0.202
#> 128  3 5.46  0.193
#> 129  3 5.52  0.185
#> 130  3 5.57  0.176
#> 131  3 5.62  0.168
#> 132  4 5.67  0.161
#> 133  4 5.74  0.153
#> 134  3 5.79  0.146
#> 135  3 5.83  0.140
#> 136  3 5.87  0.134
#> 137  4 5.90  0.128
#> 138  4 5.95  0.122
#> 139  3 5.99  0.116
#> 140  4 6.03  0.111
#> 141  4 6.07  0.106
#> 142  4 6.11  0.101
#> 143  4 6.15  0.097
#> 144  5 6.18  0.092
#> 145  5 6.24  0.088
#> 146  5 6.30  0.084
#> 147  6 6.35  0.080
#> 148  6 6.39  0.077
#> 149  6 6.43  0.073
#> 150  7 6.49  0.070
#> 151  7 6.54  0.067
#> 152  7 6.59  0.064
#> 153  7 6.64  0.061
#> 154  7 6.68  0.058
#> 155  6 6.72  0.055
#> 156  6 6.75  0.053
#> 157  7 6.79  0.051
#> 158  7 6.82  0.048
#> 159  7 6.85  0.046
#> 160  7 6.88  0.044
#> 161  8 6.91  0.042
#> 162  9 6.94  0.040
#> 163 11 6.96  0.038
#> 164 12 6.99  0.037
#> 165 13 7.02  0.035
#> 166 14 7.05  0.033
#> 167 14 7.08  0.032
#> 168 14 7.10  0.030
#> 169 15 7.13  0.029
#> 170 15 7.15  0.028
#> 171 17 7.17  0.026
#> 172 18 7.19  0.025
#> 173 20 7.21  0.024
#> 174 20 7.23  0.023
#> 175 20 7.26  0.022
#> 176 20 7.28  0.021
#> 177 21 7.30  0.020
#> 178 20 7.32  0.019
#> 179 22 7.34  0.018
#> 180 23 7.36  0.017
#> 181 23 7.38  0.017
#> 182 23 7.39  0.016
#> 183 24 7.41  0.015
#> 184 23 7.43  0.014
#> 185 23 7.44  0.014
#> 186 23 7.45  0.013
#> 187 24 7.47  0.013
#> 188 24 7.48  0.012
#> 189 26 7.49  0.011
#> 190 27 7.50  0.011
#> 191 27 7.52  0.010
#> 192 27 7.53  0.010
#> 193 28 7.54  0.010
#> 194 30 7.55  0.009
#> 195 32 7.56  0.009
#> 196 33 7.58  0.008
#> 197 33 7.59  0.008
#> 198 33 7.60  0.008
#> 199 33 7.61  0.007
#> 200 33 7.62  0.007
max_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
max_t1$performance
#> # A tibble: 3 × 33
#>   model threshold     thr_value n_presences n_absences TPR_mean  TPR_sd TNR_mean
#>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>   <dbl>    <dbl>
#> 1 max   equal_sens_s…     0.606         700        700    0.664 0.00676    0.664
#> 2 max   max_sens_spec     0.453         700        700    0.864 0.0240     0.527
#> 3 max   max_sorensen      0.442         700        700    0.924 0.0460     0.453
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
max_t1$performance_part
#> # A tibble: 9 × 21
#>   replica partition model threshold thr_value n_presences n_absences   TPR   TNR
#>   <chr>   <chr>     <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#> 1 1       1         max   max_sore…     0.368         234        234 0.953 0.397
#> 2 1       1         max   max_sens…     0.484         234        234 0.838 0.543
#> 3 1       1         max   equal_se…     0.595         234        234 0.667 0.667
#> 4 1       2         max   max_sore…     0.455         233        233 0.948 0.429
#> 5 1       2         max   max_sens…     0.523         233        233 0.884 0.506
#> 6 1       2         max   equal_se…     0.611         233        233 0.670 0.670
#> 7 1       3         max   max_sore…     0.503         233        233 0.871 0.532
#> 8 1       3         max   max_sens…     0.503         233        233 0.871 0.532
#> 9 1       3         max   equal_se…     0.606         233        233 0.657 0.657
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab   pred
#>    <chr>  <chr>      <chr> <dbl>  <dbl>
#>  1 2      .part      1         0 0.320 
#>  2 11     .part      1         0 0.340 
#>  3 12     .part      1         0 0.467 
#>  4 16     .part      1         0 0.483 
#>  5 18     .part      1         0 0.0442
#>  6 20     .part      1         0 0.115 
#>  7 23     .part      1         0 0.577 
#>  8 24     .part      1         0 0.804 
#>  9 25     .part      1         0 0.748 
#> 10 27     .part      1         0 0.665 
#> # ℹ 1,390 more rows
# }
```
