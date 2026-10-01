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
#> 1    0 0.00 1657.00
#> 2    0 0.00 1582.00
#> 3    0 0.00 1510.00
#> 4    0 0.00 1442.00
#> 5    0 0.00 1377.00
#> 6    0 0.00 1315.00
#> 7    0 0.00 1255.00
#> 8    0 0.00 1198.00
#> 9    0 0.00 1144.00
#> 10   0 0.00 1092.00
#> 11   0 0.00 1043.00
#> 12   0 0.00  995.80
#> 13   0 0.00  950.80
#> 14   0 0.00  907.80
#> 15   0 0.00  866.70
#> 16   0 0.00  827.50
#> 17   0 0.00  790.10
#> 18   0 0.00  754.40
#> 19   0 0.00  720.20
#> 20   0 0.00  687.70
#> 21   0 0.00  656.60
#> 22   0 0.00  626.90
#> 23   0 0.00  598.50
#> 24   0 0.00  571.40
#> 25   0 0.00  545.60
#> 26   0 0.00  520.90
#> 27   0 0.00  497.40
#> 28   0 0.00  474.90
#> 29   0 0.00  453.40
#> 30   0 0.00  432.90
#> 31   0 0.00  413.30
#> 32   0 0.00  394.60
#> 33   0 0.00  376.80
#> 34   0 0.00  359.70
#> 35   0 0.00  343.50
#> 36   0 0.00  327.90
#> 37   0 0.00  313.10
#> 38   0 0.00  298.90
#> 39   0 0.00  285.40
#> 40   0 0.00  272.50
#> 41   0 0.00  260.20
#> 42   0 0.00  248.40
#> 43   0 0.00  237.20
#> 44   0 0.00  226.40
#> 45   0 0.00  216.20
#> 46   0 0.00  206.40
#> 47   0 0.00  197.10
#> 48   0 0.00  188.20
#> 49   0 0.00  179.70
#> 50   0 0.00  171.50
#> 51   0 0.00  163.80
#> 52   0 0.00  156.40
#> 53   0 0.00  149.30
#> 54   0 0.00  142.50
#> 55   0 0.00  136.10
#> 56   0 0.00  129.90
#> 57   0 0.00  124.10
#> 58   0 0.00  118.50
#> 59   0 0.00  113.10
#> 60   0 0.00  108.00
#> 61   0 0.00  103.10
#> 62   0 0.00   98.44
#> 63   0 0.00   93.98
#> 64   0 0.00   89.73
#> 65   0 0.00   85.68
#> 66   0 0.00   81.80
#> 67   0 0.00   78.10
#> 68   0 0.00   74.57
#> 69   0 0.00   71.20
#> 70   0 0.00   67.98
#> 71   0 0.00   64.90
#> 72   0 0.00   61.97
#> 73   0 0.00   59.16
#> 74   0 0.00   56.49
#> 75   0 0.00   53.93
#> 76   0 0.00   51.49
#> 77   0 0.00   49.16
#> 78   0 0.00   46.94
#> 79   0 0.00   44.82
#> 80   0 0.00   42.79
#> 81   0 0.00   40.86
#> 82   0 0.00   39.01
#> 83   0 0.00   37.24
#> 84   0 0.00   35.56
#> 85   0 0.00   33.95
#> 86   0 0.00   32.41
#> 87   0 0.00   30.95
#> 88   0 0.00   29.55
#> 89   0 0.00   28.21
#> 90   0 0.00   26.94
#> 91   0 0.00   25.72
#> 92   0 0.00   24.56
#> 93   0 0.00   23.44
#> 94   0 0.00   22.38
#> 95   0 0.00   21.37
#> 96   0 0.00   20.41
#> 97   0 0.00   19.48
#> 98   0 0.00   18.60
#> 99   0 0.00   17.76
#> 100  0 0.00   16.96
#> 101  0 0.00   16.19
#> 102  0 0.00   15.46
#> 103  0 0.00   14.76
#> 104  0 0.00   14.09
#> 105  0 0.00   13.45
#> 106  0 0.00   12.85
#> 107  1 0.05   12.26
#> 108  1 0.11   11.71
#> 109  1 0.17   11.18
#> 110  1 0.22   10.67
#> 111  1 0.27   10.19
#> 112  1 0.31    9.73
#> 113  1 0.35    9.29
#> 114  1 0.38    8.87
#> 115  1 0.41    8.47
#> 116  1 0.44    8.09
#> 117  1 0.47    7.72
#> 118  1 0.49    7.37
#> 119  2 0.53    7.04
#> 120  2 0.61    6.72
#> 121  2 0.69    6.42
#> 122  2 0.77    6.12
#> 123  2 0.83    5.85
#> 124  2 0.89    5.58
#> 125  2 0.94    5.33
#> 126  2 0.99    5.09
#> 127  2 1.04    4.86
#> 128  2 1.08    4.64
#> 129  2 1.12    4.43
#> 130  2 1.16    4.23
#> 131  2 1.19    4.04
#> 132  3 1.23    3.86
#> 133  3 1.27    3.68
#> 134  3 1.31    3.52
#> 135  3 1.34    3.36
#> 136  3 1.38    3.20
#> 137  3 1.41    3.06
#> 138  4 1.44    2.92
#> 139  4 1.48    2.79
#> 140  4 1.51    2.66
#> 141  5 1.54    2.54
#> 142  5 1.58    2.43
#> 143  5 1.61    2.32
#> 144  5 1.65    2.21
#> 145  5 1.67    2.11
#> 146  5 1.70    2.02
#> 147  5 1.73    1.93
#> 148  5 1.75    1.84
#> 149  5 1.77    1.76
#> 150  5 1.79    1.68
#> 151  5 1.81    1.60
#> 152  6 1.83    1.53
#> 153  7 1.85    1.46
#> 154  7 1.88    1.39
#> 155  7 1.90    1.33
#> 156  7 1.92    1.27
#> 157  7 1.93    1.21
#> 158  7 1.95    1.16
#> 159  8 1.98    1.10
#> 160  8 2.04    1.05
#> 161  8 2.11    1.01
#> 162  8 2.16    0.96
#> 163  9 2.21    0.92
#> 164 10 2.26    0.88
#> 165 11 2.31    0.84
#> 166 11 2.36    0.80
#> 167 11 2.40    0.76
#> 168 11 2.44    0.73
#> 169 12 2.48    0.70
#> 170 13 2.52    0.66
#> 171 13 2.55    0.63
#> 172 13 2.58    0.61
#> 173 14 2.61    0.58
#> 174 14 2.64    0.55
#> 175 14 2.66    0.53
#> 176 14 2.68    0.50
#> 177 14 2.70    0.48
#> 178 15 2.72    0.46
#> 179 15 2.74    0.44
#> 180 15 2.76    0.42
#> 181 15 2.77    0.40
#> 182 15 2.79    0.38
#> 183 15 2.80    0.36
#> 184 16 2.81    0.35
#> 185 16 2.82    0.33
#> 186 16 2.84    0.32
#> 187 16 2.85    0.30
#> 188 16 2.86    0.29
#> 189 16 2.87    0.28
#> 190 16 2.88    0.26
#> 191 17 2.89    0.25
#> 192 17 2.90    0.24
#> 193 17 2.91    0.23
#> 194 17 2.92    0.22
#> 195 17 2.93    0.21
#> 196 17 2.93    0.20
#> 197 17 2.94    0.19
#> 198 16 2.95    0.18
#> 199 16 2.95    0.17
#> 200 16 2.96    0.17
max_t1$predictors
#> # A tibble: 1 × 5
#>   c1    c2    c3    c4    f       
#>   <chr> <chr> <chr> <chr> <chr>   
#> 1 aet   pH    awc   depth landform
max_t1$performance
#> # A tibble: 1 × 35
#>   regmult classes model threshold     thr_value n_presences n_absences TPR_mean
#>     <dbl> <fct>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>
#> 1       1 lq      max   max_sens_spec     0.522         700        700    0.844
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
#> 1 1       1         max   max_sens…     0.503         234        234 0.872 0.521
#> 2 1       2         max   max_sens…     0.575         233        233 0.773 0.579
#> 3 1       3         max   max_sens…     0.490         233        233 0.888 0.536
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <dbl> <dbl>
#>  1 1      .part      1         0 0.563
#>  2 3      .part      1         0 0.671
#>  3 5      .part      1         0 0.666
#>  4 7      .part      1         0 0.424
#>  5 8      .part      1         0 0.747
#>  6 13     .part      1         0 0.655
#>  7 15     .part      1         0 0.689
#>  8 16     .part      1         0 0.487
#>  9 22     .part      1         0 0.326
#> 10 23     .part      1         0 0.454
#> # ℹ 1,390 more rows
# }
```
