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
# Note that the partition method, number of folds or replications must
# be the same for presence-absence and background points datasets
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
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
  method = c(method = "kfold", folds = 5)
)
backg2
#> # A tibble: 5,000 × 14
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
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5
length(max_t1)
#> [1] 5

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
#> 7    0  0.00 16.1900
#> 8    0  0.00 15.4600
#> 9    0  0.00 14.7600
#> 10   0  0.00 14.0900
#> 11   0  0.00 13.4500
#> 12   0  0.00 12.8500
#> 13   0  0.00 12.2600
#> 14   0  0.00 11.7100
#> 15   0  0.00 11.1800
#> 16   0  0.00 10.6700
#> 17   0  0.00 10.1900
#> 18   0  0.00  9.7310
#> 19   0  0.00  9.2900
#> 20   0  0.00  8.8700
#> 21   0  0.00  8.4690
#> 22   0  0.00  8.0860
#> 23   0  0.00  7.7200
#> 24   0  0.00  7.3710
#> 25   0  0.00  7.0380
#> 26   0  0.00  6.7190
#> 27   0  0.00  6.4160
#> 28   0  0.00  6.1250
#> 29   0  0.00  5.8480
#> 30   0  0.00  5.5840
#> 31   0  0.00  5.3310
#> 32   0  0.00  5.0900
#> 33   0  0.00  4.8600
#> 34   0  0.00  4.6400
#> 35   0  0.00  4.4300
#> 36   0  0.00  4.2300
#> 37   0  0.00  4.0390
#> 38   0  0.00  3.8560
#> 39   0  0.00  3.6820
#> 40   0  0.00  3.5150
#> 41   0  0.00  3.3560
#> 42   0  0.00  3.2040
#> 43   0  0.00  3.0590
#> 44   0  0.00  2.9210
#> 45   0  0.00  2.7890
#> 46   0  0.00  2.6630
#> 47   0  0.00  2.5420
#> 48   0  0.00  2.4270
#> 49   0  0.00  2.3180
#> 50   0  0.00  2.2130
#> 51   0  0.00  2.1130
#> 52   0  0.00  2.0170
#> 53   0  0.00  1.9260
#> 54   0  0.00  1.8390
#> 55   0  0.00  1.7560
#> 56   0  0.00  1.6760
#> 57   0  0.00  1.6000
#> 58   0  0.00  1.5280
#> 59   0  0.00  1.4590
#> 60   0  0.00  1.3930
#> 61   0  0.00  1.3300
#> 62   0  0.00  1.2700
#> 63   0  0.00  1.2120
#> 64   0  0.00  1.1570
#> 65   0  0.00  1.1050
#> 66   0  0.00  1.0550
#> 67   0  0.00  1.0070
#> 68   0  0.00  0.9619
#> 69   0  0.00  0.9184
#> 70   0  0.00  0.8768
#> 71   0  0.00  0.8372
#> 72   0  0.00  0.7993
#> 73   0  0.00  0.7632
#> 74   0  0.00  0.7286
#> 75   1  0.21  0.6957
#> 76   1  0.65  0.6642
#> 77   1  1.06  0.6342
#> 78   1  1.44  0.6055
#> 79   1  1.80  0.5781
#> 80   1  2.13  0.5520
#> 81   1  2.45  0.5270
#> 82   1  2.74  0.5032
#> 83   1  3.01  0.4804
#> 84   1  3.27  0.4587
#> 85   1  3.52  0.4379
#> 86   1  3.75  0.4181
#> 87   1  3.96  0.3992
#> 88   1  4.17  0.3812
#> 89   1  4.36  0.3639
#> 90   1  4.54  0.3475
#> 91   1  4.72  0.3317
#> 92   1  4.88  0.3167
#> 93   1  5.04  0.3024
#> 94   1  5.19  0.2887
#> 95   1  5.33  0.2757
#> 96   1  5.46  0.2632
#> 97   1  5.59  0.2513
#> 98   1  5.71  0.2399
#> 99   2  5.88  0.2291
#> 100  2  6.07  0.2187
#> 101  2  6.25  0.2088
#> 102  2  6.42  0.1994
#> 103  2  6.58  0.1904
#> 104  3  6.75  0.1818
#> 105  3  6.93  0.1735
#> 106  3  7.11  0.1657
#> 107  3  7.27  0.1582
#> 108  3  7.42  0.1510
#> 109  3  7.56  0.1442
#> 110  3  7.69  0.1377
#> 111  3  7.81  0.1315
#> 112  3  7.92  0.1255
#> 113  3  8.03  0.1198
#> 114  3  8.12  0.1144
#> 115  3  8.22  0.1092
#> 116  2  8.30  0.1043
#> 117  2  8.38  0.0996
#> 118  2  8.45  0.0951
#> 119  2  8.51  0.0908
#> 120  4  8.58  0.0867
#> 121  4  8.71  0.0828
#> 122  5  8.84  0.0790
#> 123  5  9.00  0.0754
#> 124  5  9.14  0.0720
#> 125  5  9.28  0.0688
#> 126  6  9.44  0.0657
#> 127  6  9.59  0.0627
#> 128  6  9.74  0.0598
#> 129  6  9.88  0.0572
#> 130  6 10.01  0.0546
#> 131  7 10.13  0.0521
#> 132  6 10.22  0.0497
#> 133  6 10.31  0.0475
#> 134  6 10.39  0.0453
#> 135  7 10.48  0.0433
#> 136  7 10.60  0.0413
#> 137  6 10.71  0.0395
#> 138  6 10.81  0.0377
#> 139  7 10.90  0.0360
#> 140  7 10.98  0.0344
#> 141  7 11.06  0.0328
#> 142  7 11.13  0.0313
#> 143  7 11.20  0.0299
#> 144  8 11.26  0.0285
#> 145  9 11.33  0.0272
#> 146  9 11.40  0.0260
#> 147  9 11.46  0.0248
#> 148  9 11.52  0.0237
#> 149 10 11.58  0.0226
#> 150 10 11.64  0.0216
#> 151 11 11.70  0.0206
#> 152 11 11.77  0.0197
#> 153 11 11.82  0.0188
#> 154 11 11.88  0.0180
#> 155 12 11.93  0.0171
#> 156 15 11.99  0.0164
#> 157 16 12.05  0.0156
#> 158 18 12.12  0.0149
#> 159 20 12.18  0.0143
#> 160 20 12.25  0.0136
#> 161 20 12.31  0.0130
#> 162 21 12.36  0.0124
#> 163 21 12.42  0.0118
#> 164 21 12.46  0.0113
#> 165 21 12.51  0.0108
#> 166 21 12.55  0.0103
#> 167 20 12.59  0.0098
#> 168 20 12.63  0.0094
#> 169 21 12.67  0.0090
#> 170 21 12.70  0.0086
#> 171 22 12.73  0.0082
#> 172 22 12.76  0.0078
#> 173 24 12.80  0.0075
#> 174 24 12.83  0.0071
#> 175 22 12.86  0.0068
#> 176 24 12.89  0.0065
#> 177 24 12.92  0.0062
#> 178 24 12.95  0.0059
#> 179 25 12.97  0.0056
#> 180 24 13.00  0.0054
#> 181 25 13.02  0.0051
#> 182 26 13.05  0.0049
#> 183 27 13.08  0.0047
#> 184 27 13.11  0.0045
#> 185 27 13.14  0.0043
#> 186 26 13.15  0.0041
#> 187 28 13.17  0.0039
#> 188 27 13.20  0.0037
#> 189 27 13.23  0.0036
#> 190 30 13.25  0.0034
#> 191 32 13.28  0.0032
#> 192 33 13.30  0.0031
#> 193 34 13.33  0.0030
#> 194 34 13.35  0.0028
#> 195 34 13.37  0.0027
#> 196 35 13.40  0.0026
#> 197 36 13.43  0.0025
#> 198 37 13.46  0.0023
#> 199 40 13.51  0.0022
#> 200 43 13.55  0.0021
max_t1$predictors
#> # A tibble: 1 × 6
#>   c1    c2      c3    c4    c5    f       
#>   <chr> <chr>   <chr> <chr> <chr> <chr>   
#> 1 aet   ppt_jja pH    awc   depth landform
max_t1$performance
#> # A tibble: 3 × 33
#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean
#>   <chr> <chr>              <dbl>       <int>      <int>    <dbl>  <dbl>    <dbl>
#> 1 max   equal_sens_sp…     0.573         700        700    0.667 0.0393    0.667
#> 2 max   max_sens_spec      0.416         700        700    0.874 0.0690    0.55 
#> 3 max   max_sorensen       0.336         700        700    0.934 0.0287    0.477
#> # ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,
#> #   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,
#> #   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,
#> #   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,
#> #   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,
#> #   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,
#> #   IMAE_mean <dbl>, IMAE_sd <dbl>
max_t1$performance_part
#> # A tibble: 15 × 21
#>    replica partition model threshold      thr_value n_presences n_absences   TPR
#>    <chr>   <chr>     <chr> <chr>              <dbl>       <int>      <int> <dbl>
#>  1 1       5         max   max_sorensen       0.322         140        140 0.95 
#>  2 1       5         max   max_sens_spec      0.507         140        140 0.786
#>  3 1       5         max   equal_sens_sp…     0.565         140        140 0.714
#>  4 1       2         max   max_sorensen       0.355         140        140 0.9  
#>  5 1       2         max   max_sens_spec      0.355         140        140 0.9  
#>  6 1       2         max   equal_sens_sp…     0.541         140        140 0.664
#>  7 1       3         max   max_sorensen       0.316         140        140 0.95 
#>  8 1       3         max   max_sens_spec      0.487         140        140 0.829
#>  9 1       3         max   equal_sens_sp…     0.620         140        140 0.607
#> 10 1       4         max   max_sorensen       0.316         140        140 0.907
#> 11 1       4         max   max_sens_spec      0.322         140        140 0.893
#> 12 1       4         max   equal_sens_sp…     0.561         140        140 0.664
#> 13 1       5         max   max_sorensen       0.336         140        140 0.964
#> 14 1       5         max   max_sens_spec      0.336         140        140 0.964
#> 15 1       5         max   equal_sens_sp…     0.562         140        140 0.686
#> # ℹ 13 more variables: TNR <dbl>, W_TPR_TNR <dbl>, SORENSEN <dbl>,
#> #   JACCARD <dbl>, FPB <dbl>, OR <dbl>, TSS <dbl>, KAPPA <dbl>, MCC <dbl>,
#> #   AUC <dbl>, BOYCE <dbl>, CRPS <dbl>, IMAE <dbl>
max_t1$data_ens
#> # A tibble: 1,400 × 5
#>    rnames replicates part  pr_ab   pred
#>    <chr>  <chr>      <chr> <dbl>  <dbl>
#>  1 16     .part      1         0 0.261 
#>  2 17     .part      1         0 0.481 
#>  3 23     .part      1         0 0.523 
#>  4 30     .part      1         0 0.0174
#>  5 31     .part      1         0 0.404 
#>  6 34     .part      1         0 0.309 
#>  7 38     .part      1         0 0.189 
#>  8 42     .part      1         0 0.248 
#>  9 43     .part      1         0 0.901 
#> 10 48     .part      1         0 0.153 
#> # ℹ 1,390 more rows
# }
```
