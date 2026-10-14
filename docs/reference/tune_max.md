# Fit and validate Maximum Entropy models with exploration of hyper-parameters that optimize performance

Fit and validate Maximum Entropy models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_max(
  data,
  response,
  predictors,
  predictors_f = NULL,
  background = NULL,
  partition,
  grid = NULL,
  thr = NULL,
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 1
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

- background:

  data.frame. Database with response variable column only containing 0
  values, and predictors variables. All column names must be consistent
  with data

- partition:

  character. Column name with training and validation partition groups.

- grid:

  data.frame. A data frame object with algorithm hyper-parameters values
  to be tested. It is recommended to generate this data.frame with the
  grid() function. Hyper-parameters needed for tuning are 'regmult' and
  'classes' (any combination of following letters l -linear-, q
  -quadratic-, h -hinge-, p -product-, and t -threshold-).

- thr:

  character. Threshold used to get binary suitability values (i.e.
  0,1)., needed for threshold-dependent performance metrics. More than
  one threshold type can be used. It is necessary to provide a vector
  for this argument. The following threshold types are available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which sensitivity and specificity are
    equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which \# FPB (F-measure on
    presence-background data) is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If a sensitivity value is not
    specified, a default of 0.9 will be used.

  If more than one threshold type is used, concatenate them, e.g.,
  thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold is specified.

- metric:

  character. Performance metric used for selecting the best combination
  of hyper -parameter values. One of the following metrics can be used:
  SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as
  default.

- clamp:

  logical. If TRUE, predictors and features are restricted to the range
  seen during model training.

- pred_type:

  character. Type of response required available "link", "exponential",
  "cloglog" and "logistic". Default "cloglog"

- n_cores:

  numeric. Number of cores use for parallelization. Default 1

## Value

A list object with:

- model: A "maxnet" class object from maxnet package. This object can be
  used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Hyper-parameters values and performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md))
  for the best hyper-parameters combination.

- performance_part: Performance metric for each replica and partition
  (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)).

- hyper_performance: Performance metrics (see
  [`sdm_eval`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md))
  for each combination of the hyper-parameters.

- data_ens: Predicted suitability for each test partition based on the
  best model. This database is used in
  [`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)

## Details

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

[`tune_gbm`](https://sjevelazco.github.io/flexsdm/reference/tune_gbm.md),
[`tune_net`](https://sjevelazco.github.io/flexsdm/reference/tune_net.md),
[`tune_raf`](https://sjevelazco.github.io/flexsdm/reference/tune_raf.md),
and
[`tune_svm`](https://sjevelazco.github.io/flexsdm/reference/tune_svm.md).

## Examples

``` r
# Toy data: a small subset of 'abies' to keep this example fast
data("abies")
abies_toy <- abies[c(1:30, 701:730), ]
abies_toy <- part_random(
  data = abies_toy,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

# Use the same records as background points (toy example only)
backg_toy <- abies_toy
backg_toy$pr_ab <- 0

tune_toy <- tune_max(
  data = abies_toy,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmin"),
  background = backg_toy,
  partition = ".part",
  grid = expand.grid(regmult = c(0.5, 1), classes = c("l", "lq")),
  thr = "max_sens_spec",
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 1
)
#> Tuning model...
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Fitting best model
#> Formula used for model fitting:
#> ~aet + cwd + tmin + I(aet^2) + I(cwd^2) + I(tmin^2) - 1
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
tune_toy$performance
#> # A tibble: 1 × 35
#>   regmult classes model threshold     thr_value n_presences n_absences TPR_mean
#>     <dbl> <fct>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>
#> 1     0.5 lq      max   max_sens_spec     0.397          30         30    0.933
#> # ℹ 27 more variables: TPR_sd <dbl>, TNR_mean <dbl>, TNR_sd <dbl>,
#> #   W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>,
#> #   SORENSEN_sd <dbl>, JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>,
#> #   FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>,
#> #   KAPPA_mean <dbl>, KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>,
#> #   AUC_mean <dbl>, AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>,
#> #   CRPS_mean <dbl>, CRPS_sd <dbl>, IMAE_mean <dbl>, IMAE_sd <dbl>

# \donttest{
data("abies")
data("backg")
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
#> # A tibble: 5,000 × 13
#>    pr_ab        x        y   aet   cwd  tmin ppt_djf ppt_jja    pH     awc depth
#>    <dbl>    <dbl>    <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl>   <dbl> <dbl>
#>  1     0  160779. -449968.  280. 1137. 13.5     71.3    1.19 0     0         0  
#>  2     0   36849.   24152.  260.  382. -3.17   171.    17.5  0.212 0.00347 201  
#>  3     0 -240171.   90032.  400.  700.  8.68   285.     5.02 5.72  0.0804   50.1
#>  4     0 -152421. -143518.  367.  843.  9.01    72.0    1.20 7.54  0.170   154. 
#>  5     0 -193191.   24152.  397.  842.  8.97   125.     1.98 6.20  0.131   122. 
#>  6     0 -277971.  223682.  385.  637.  4.93   226.     8.16 5.81  0.0512   56.2
#>  7     0 -313341.  270122.  582.  406.  6.29   334.    18.4  5.80  0.168   201  
#>  8     0   54399.  -15538.  346.  195. -5.22   142.    13.0  5.60  0.120   201  
#>  9     0  282549. -582268.  285. 1097. 11.4     56.3    1.39 6.57  0.135    68.4
#> 10     0  104079. -178618.  385.  871.  6.76   147.     7.80 6.10  0.0300   41  
#> # ℹ 4,990 more rows
#> # ℹ 2 more variables: percent_clay <dbl>, landform <fct>

# Using k-fold partition method
# Remember that the partition method, number of folds or replications must
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

set.seed(1)
backg <- dplyr::sample_n(backg, size = 500, replace = FALSE)
backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
backg
#> # A tibble: 500 × 13
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
#> # ℹ 490 more rows
#> # ℹ 2 more variables: percent_clay <dbl>, landform <fct>


gridtest <-
  expand.grid(
    regmult = c(0.1, 1),
    classes = c("l", "lq")
  )

max_t1 <- tune_max(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  background = backg2,
  grid = gridtest,
  thr = "max_sens_spec",
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 2 # activate two cores to speed up this process
)
#> Tuning model...
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3
#> Fitting best model
#> Formula used for model fitting:
#> ~aet + pH + awc + depth + I(aet^2) + I(pH^2) + I(awc^2) + I(depth^2) + categorical(landform) - 1
#> Replica number: 1/1
#> Partition number: 1/3
#> Partition number: 2/3
#> Partition number: 3/3

length(max_t1)
#> [1] 6
max_t1$model
#> 
#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) 
#> 
#>     Df %Dev  Lambda
#> 1    0 0.00 165.700
#> 2    0 0.00 158.200
#> 3    0 0.00 151.000
#> 4    0 0.00 144.200
#> 5    0 0.00 137.700
#> 6    0 0.00 131.500
#> 7    0 0.00 125.500
#> 8    0 0.00 119.800
#> 9    0 0.00 114.400
#> 10   0 0.00 109.200
#> 11   0 0.00 104.300
#> 12   0 0.00  99.580
#> 13   0 0.00  95.080
#> 14   0 0.00  90.780
#> 15   0 0.00  86.670
#> 16   0 0.00  82.750
#> 17   0 0.00  79.010
#> 18   0 0.00  75.440
#> 19   0 0.00  72.020
#> 20   0 0.00  68.770
#> 21   0 0.00  65.660
#> 22   0 0.00  62.690
#> 23   0 0.00  59.850
#> 24   0 0.00  57.140
#> 25   0 0.00  54.560
#> 26   0 0.00  52.090
#> 27   0 0.00  49.740
#> 28   0 0.00  47.490
#> 29   0 0.00  45.340
#> 30   0 0.00  43.290
#> 31   0 0.00  41.330
#> 32   0 0.00  39.460
#> 33   0 0.00  37.680
#> 34   0 0.00  35.970
#> 35   0 0.00  34.350
#> 36   0 0.00  32.790
#> 37   0 0.00  31.310
#> 38   0 0.00  29.890
#> 39   0 0.00  28.540
#> 40   0 0.00  27.250
#> 41   0 0.00  26.020
#> 42   0 0.00  24.840
#> 43   0 0.00  23.720
#> 44   0 0.00  22.640
#> 45   0 0.00  21.620
#> 46   0 0.00  20.640
#> 47   0 0.00  19.710
#> 48   0 0.00  18.820
#> 49   0 0.00  17.970
#> 50   0 0.00  17.150
#> 51   0 0.00  16.380
#> 52   0 0.00  15.640
#> 53   0 0.00  14.930
#> 54   0 0.00  14.250
#> 55   0 0.00  13.610
#> 56   0 0.00  12.990
#> 57   1 0.03  12.410
#> 58   1 0.10  11.850
#> 59   1 0.15  11.310
#> 60   1 0.21  10.800
#> 61   1 0.26  10.310
#> 62   1 0.30   9.844
#> 63   1 0.34   9.398
#> 64   1 0.37   8.973
#> 65   1 0.41   8.568
#> 66   1 0.43   8.180
#> 67   1 0.46   7.810
#> 68   1 0.48   7.457
#> 69   2 0.51   7.120
#> 70   2 0.59   6.798
#> 71   2 0.67   6.490
#> 72   3 0.74   6.197
#> 73   2 0.82   5.916
#> 74   2 0.88   5.649
#> 75   2 0.93   5.393
#> 76   2 0.98   5.149
#> 77   2 1.03   4.916
#> 78   2 1.07   4.694
#> 79   2 1.11   4.482
#> 80   2 1.15   4.279
#> 81   2 1.18   4.086
#> 82   3 1.22   3.901
#> 83   3 1.26   3.724
#> 84   3 1.30   3.556
#> 85   3 1.33   3.395
#> 86   3 1.37   3.241
#> 87   3 1.40   3.095
#> 88   4 1.43   2.955
#> 89   4 1.47   2.821
#> 90   4 1.50   2.694
#> 91   5 1.53   2.572
#> 92   5 1.57   2.456
#> 93   5 1.61   2.344
#> 94   5 1.64   2.238
#> 95   5 1.67   2.137
#> 96   5 1.70   2.041
#> 97   5 1.72   1.948
#> 98   5 1.74   1.860
#> 99   5 1.77   1.776
#> 100  5 1.79   1.696
#> 101  5 1.81   1.619
#> 102  6 1.82   1.546
#> 103  7 1.85   1.476
#> 104  7 1.87   1.409
#> 105  7 1.89   1.345
#> 106  7 1.91   1.285
#> 107  7 1.93   1.226
#> 108  7 1.95   1.171
#> 109  8 1.97   1.118
#> 110  8 2.03   1.067
#> 111  8 2.09   1.019
#> 112  8 2.15   0.973
#> 113  8 2.20   0.929
#> 114 10 2.25   0.887
#> 115 11 2.30   0.847
#> 116 11 2.35   0.809
#> 117 11 2.39   0.772
#> 118 11 2.43   0.737
#> 119 12 2.47   0.704
#> 120 13 2.51   0.672
#> 121 13 2.54   0.641
#> 122 13 2.57   0.613
#> 123 14 2.60   0.585
#> 124 14 2.63   0.558
#> 125 14 2.65   0.533
#> 126 14 2.68   0.509
#> 127 14 2.70   0.486
#> 128 15 2.72   0.464
#> 129 15 2.74   0.443
#> 130 15 2.75   0.423
#> 131 15 2.77   0.404
#> 132 15 2.78   0.386
#> 133 15 2.80   0.368
#> 134 16 2.81   0.352
#> 135 16 2.82   0.336
#> 136 16 2.83   0.320
#> 137 16 2.85   0.306
#> 138 16 2.86   0.292
#> 139 16 2.87   0.279
#> 140 16 2.88   0.266
#> 141 16 2.89   0.254
#> 142 17 2.90   0.243
#> 143 17 2.91   0.232
#> 144 17 2.92   0.221
#> 145 17 2.92   0.211
#> 146 17 2.93   0.202
#> 147 17 2.94   0.193
#> 148 16 2.95   0.184
#> 149 16 2.95   0.176
#> 150 16 2.96   0.168
#> 151 16 2.96   0.160
#> 152 16 2.97   0.153
#> 153 17 2.97   0.146
#> 154 17 2.97   0.139
#> 155 17 2.98   0.133
#> 156 18 2.98   0.127
#> 157 18 2.98   0.121
#> 158 18 2.99   0.116
#> 159 18 2.99   0.110
#> 160 18 2.99   0.106
#> 161 18 2.99   0.101
#> 162 18 2.99   0.096
#> 163 19 3.00   0.092
#> 164 19 3.00   0.088
#> 165 20 3.00   0.084
#> 166 20 3.00   0.080
#> 167 20 3.00   0.076
#> 168 20 3.01   0.073
#> 169 20 3.01   0.070
#> 170 20 3.01   0.066
#> 171 20 3.01   0.063
#> 172 19 3.01   0.061
#> 173 19 3.01   0.058
#> 174 19 3.01   0.055
#> 175 19 3.01   0.053
#> 176 19 3.01   0.050
#> 177 19 3.01   0.048
#> 178 19 3.01   0.046
#> 179 19 3.01   0.044
#> 180 19 3.01   0.042
#> 181 19 3.01   0.040
#> 182 19 3.02   0.038
#> 183 19 3.02   0.036
#> 184 19 3.02   0.035
#> 185 19 3.02   0.033
#> 186 19 3.02   0.032
#> 187 19 3.02   0.030
#> 188 19 3.02   0.029
#> 189 19 3.02   0.028
#> 190 19 3.02   0.026
#> 191 19 3.02   0.025
#> 192 19 3.02   0.024
#> 193 19 3.02   0.023
#> 194 19 3.02   0.022
#> 195 19 3.02   0.021
#> 196 20 3.02   0.020
#> 197 20 3.02   0.019
#> 198 20 3.02   0.018
#> 199 20 3.02   0.017
#> 200 20 3.02   0.017
max_t1$predictors
#> # A tibble: 1 × 5
#>   c1    c2    c3    c4    f       
#>   <chr> <chr> <chr> <chr> <chr>   
#> 1 aet   pH    awc   depth landform
max_t1$performance
#> # A tibble: 1 × 35
#>   regmult classes model threshold     thr_value n_presences n_absences TPR_mean
#>     <dbl> <fct>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>
#> 1     0.1 lq      max   max_sens_spec     0.516         700        700    0.830
#> # ℹ 27 more variables: TPR_sd <dbl>, TNR_mean <dbl>, TNR_sd <dbl>,
#> #   W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>,
#> #   SORENSEN_sd <dbl>, JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>,
#> #   FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>,
#> #   KAPPA_mean <dbl>, KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>,
#> #   AUC_mean <dbl>, AUC_sd <dbl>, BOYCE_mean <dbl>, BOYCE_sd <dbl>,
#> #   CRPS_mean <dbl>, CRPS_sd <dbl>, IMAE_mean <dbl>, IMAE_sd <dbl>
max_t1$performance_part
#> # A tibble: 3 × 21
#>   replica partition model threshold thr_value n_presences n_absences   TPR   TNR
#>   <chr>   <chr>     <chr> <chr>         <dbl>       <int>      <int> <dbl> <dbl>
#> 1 1       1         max   max_sens…     0.520         234        234 0.821 0.534
#> 2 1       2         max   max_sens…     0.530         233        233 0.854 0.605
#> 3 1       3         max   max_sens…     0.508         233        233 0.815 0.515
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <dbl> <dbl>
#>  1 3      .part      1         0 0.535
#>  2 4      .part      1         0 0.387
#>  3 5      .part      1         0 0.530
#>  4 7      .part      1         0 0.433
#>  5 9      .part      1         0 0.497
#>  6 12     .part      1         0 0.301
#>  7 25     .part      1         0 0.720
#>  8 26     .part      1         0 0.716
#>  9 28     .part      1         0 0.649
#> 10 33     .part      1         0 0.597
#> # ℹ 1,390 more rows
# }
```
