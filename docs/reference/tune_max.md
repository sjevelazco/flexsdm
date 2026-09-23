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
backg <- dplyr::sample_n(backg, size = 2000, replace = FALSE)
backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
backg
#> # A tibble: 2,000 × 13
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
#> # ℹ 1,990 more rows
#> # ℹ 2 more variables: percent_clay <dbl>, landform <fct>


gridtest <-
  expand.grid(
    regmult = seq(0.1, 3, 0.5),
    classes = c("l", "lq", "lqh")
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
#> 1    0 0.00 1991.00
#> 2    0 0.00 1901.00
#> 3    0 0.00 1815.00
#> 4    0 0.00 1733.00
#> 5    0 0.00 1655.00
#> 6    0 0.00 1580.00
#> 7    0 0.00 1509.00
#> 8    0 0.00 1440.00
#> 9    0 0.00 1375.00
#> 10   0 0.00 1313.00
#> 11   0 0.00 1254.00
#> 12   0 0.00 1197.00
#> 13   0 0.00 1143.00
#> 14   0 0.00 1091.00
#> 15   0 0.00 1042.00
#> 16   0 0.00  994.60
#> 17   0 0.00  949.60
#> 18   0 0.00  906.70
#> 19   0 0.00  865.70
#> 20   0 0.00  826.50
#> 21   0 0.00  789.10
#> 22   0 0.00  753.50
#> 23   0 0.00  719.40
#> 24   0 0.00  686.80
#> 25   0 0.00  655.80
#> 26   0 0.00  626.10
#> 27   0 0.00  597.80
#> 28   0 0.00  570.80
#> 29   0 0.00  544.90
#> 30   0 0.00  520.30
#> 31   0 0.00  496.80
#> 32   0 0.00  474.30
#> 33   0 0.00  452.80
#> 34   0 0.00  432.40
#> 35   0 0.00  412.80
#> 36   0 0.00  394.10
#> 37   0 0.00  376.30
#> 38   0 0.00  359.30
#> 39   0 0.00  343.00
#> 40   0 0.00  327.50
#> 41   0 0.00  312.70
#> 42   0 0.00  298.60
#> 43   0 0.00  285.10
#> 44   0 0.00  272.20
#> 45   0 0.00  259.90
#> 46   0 0.00  248.10
#> 47   0 0.00  236.90
#> 48   0 0.00  226.20
#> 49   0 0.00  215.90
#> 50   0 0.00  206.20
#> 51   0 0.00  196.90
#> 52   0 0.00  187.90
#> 53   0 0.00  179.40
#> 54   0 0.00  171.30
#> 55   0 0.00  163.60
#> 56   0 0.00  156.20
#> 57   0 0.00  149.10
#> 58   0 0.00  142.40
#> 59   0 0.00  135.90
#> 60   0 0.00  129.80
#> 61   0 0.00  123.90
#> 62   0 0.00  118.30
#> 63   0 0.00  113.00
#> 64   0 0.00  107.90
#> 65   0 0.00  103.00
#> 66   0 0.00   98.32
#> 67   0 0.00   93.87
#> 68   0 0.00   89.63
#> 69   0 0.00   85.57
#> 70   0 0.00   81.70
#> 71   0 0.00   78.01
#> 72   0 0.00   74.48
#> 73   0 0.00   71.11
#> 74   0 0.00   67.89
#> 75   0 0.00   64.82
#> 76   0 0.00   61.89
#> 77   0 0.00   59.09
#> 78   0 0.00   56.42
#> 79   0 0.00   53.87
#> 80   0 0.00   51.43
#> 81   0 0.00   49.11
#> 82   0 0.00   46.88
#> 83   0 0.00   44.76
#> 84   0 0.00   42.74
#> 85   0 0.00   40.81
#> 86   0 0.00   38.96
#> 87   0 0.00   37.20
#> 88   0 0.00   35.52
#> 89   0 0.00   33.91
#> 90   0 0.00   32.38
#> 91   0 0.00   30.91
#> 92   0 0.00   29.51
#> 93   0 0.00   28.18
#> 94   0 0.00   26.90
#> 95   0 0.00   25.69
#> 96   0 0.00   24.53
#> 97   0 0.00   23.42
#> 98   0 0.00   22.36
#> 99   0 0.00   21.35
#> 100  0 0.00   20.38
#> 101  0 0.00   19.46
#> 102  0 0.00   18.58
#> 103  0 0.00   17.74
#> 104  0 0.00   16.94
#> 105  0 0.00   16.17
#> 106  0 0.00   15.44
#> 107  0 0.00   14.74
#> 108  0 0.00   14.07
#> 109  0 0.00   13.44
#> 110  0 0.00   12.83
#> 111  0 0.00   12.25
#> 112  0 0.00   11.70
#> 113  0 0.00   11.17
#> 114  0 0.00   10.66
#> 115  0 0.00   10.18
#> 116  1 0.18    9.72
#> 117  1 0.37    9.28
#> 118  1 0.54    8.86
#> 119  1 0.69    8.46
#> 120  1 0.82    8.08
#> 121  1 0.94    7.71
#> 122  1 1.05    7.36
#> 123  1 1.14    7.03
#> 124  1 1.23    6.71
#> 125  1 1.30    6.41
#> 126  1 1.37    6.12
#> 127  1 1.43    5.84
#> 128  1 1.49    5.58
#> 129  1 1.54    5.32
#> 130  1 1.58    5.08
#> 131  1 1.62    4.85
#> 132  2 1.70    4.63
#> 133  3 1.86    4.42
#> 134  3 2.11    4.22
#> 135  3 2.36    4.03
#> 136  2 2.55    3.85
#> 137  2 2.65    3.68
#> 138  2 2.75    3.51
#> 139  2 2.83    3.35
#> 140  2 2.91    3.20
#> 141  3 2.99    3.06
#> 142  3 3.08    2.92
#> 143  3 3.16    2.79
#> 144  3 3.24    2.66
#> 145  3 3.31    2.54
#> 146  3 3.38    2.42
#> 147  3 3.44    2.32
#> 148  3 3.50    2.21
#> 149  3 3.55    2.11
#> 150  3 3.60    2.02
#> 151  3 3.65    1.92
#> 152  3 3.69    1.84
#> 153  4 3.74    1.75
#> 154  4 3.78    1.67
#> 155  4 3.83    1.60
#> 156  4 3.87    1.53
#> 157  4 3.91    1.46
#> 158  4 3.95    1.39
#> 159  4 3.98    1.33
#> 160  4 4.01    1.27
#> 161  4 4.04    1.21
#> 162  4 4.07    1.16
#> 163  4 4.09    1.10
#> 164  4 4.11    1.05
#> 165  4 4.14    1.01
#> 166  5 4.16    0.96
#> 167  6 4.18    0.92
#> 168  7 4.22    0.88
#> 169  7 4.25    0.84
#> 170  7 4.29    0.80
#> 171  8 4.34    0.76
#> 172  8 4.48    0.73
#> 173  8 4.61    0.69
#> 174  9 4.73    0.66
#> 175 10 4.85    0.63
#> 176 10 4.95    0.60
#> 177 10 5.05    0.58
#> 178 10 5.13    0.55
#> 179 10 5.22    0.53
#> 180 11 5.29    0.50
#> 181 11 5.36    0.48
#> 182 11 5.42    0.46
#> 183 11 5.48    0.44
#> 184 11 5.54    0.42
#> 185 12 5.59    0.40
#> 186 13 5.63    0.38
#> 187 13 5.68    0.36
#> 188 13 5.72    0.35
#> 189 13 5.75    0.33
#> 190 14 5.79    0.32
#> 191 14 5.82    0.30
#> 192 15 5.85    0.29
#> 193 15 5.88    0.28
#> 194 15 5.90    0.26
#> 195 15 5.92    0.25
#> 196 15 5.95    0.24
#> 197 15 5.97    0.23
#> 198 15 5.98    0.22
#> 199 15 6.00    0.21
#> 200 15 6.01    0.20
max_t1$predictors
#> # A tibble: 1 × 5
#>   c1    c2    c3    c4    f       
#>   <chr> <chr> <chr> <chr> <chr>   
#> 1 aet   pH    awc   depth landform
max_t1$performance
#> # A tibble: 1 × 35
#>   regmult classes model threshold     thr_value n_presences n_absences TPR_mean
#>     <dbl> <fct>   <chr> <chr>             <dbl>       <int>      <int>    <dbl>
#> 1     2.6 lq      max   max_sens_spec     0.445         700        700    0.900
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
#> 1 1       3         max   max_sens…     0.403         234        234 0.893 0.466
#> 2 1       2         max   max_sens…     0.461         233        233 0.876 0.468
#> 3 1       3         max   max_sens…     0.385         233        233 0.931 0.464
#> # ℹ 12 more variables: W_TPR_TNR <dbl>, SORENSEN <dbl>, JACCARD <dbl>,
#> #   FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>, AUC <dbl>,
#> #   BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab  pred
#>    <chr>  <chr>      <chr> <dbl> <dbl>
#>  1 2      .part      1         0 0.212
#>  2 4      .part      1         0 0.326
#>  3 5      .part      1         0 0.796
#>  4 12     .part      1         0 0.241
#>  5 17     .part      1         0 0.523
#>  6 24     .part      1         0 0.955
#>  7 25     .part      1         0 0.817
#>  8 29     .part      1         0 0.287
#>  9 30     .part      1         0 0.529
#> 10 31     .part      1         0 0.390
#> # ℹ 1,390 more rows
# }
```
