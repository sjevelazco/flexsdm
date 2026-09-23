# Conventional data partitioning methods

This function provides different conventional (randomized, non-spatial)
partitioning methods based on cross validation folds (kfold, rep_kfold,
and loocv), as well as bootstrap (boot)

## Usage

``` r
part_random(data, pr_ab, method = NULL)
```

## Arguments

- data:

  data.frame. Database with presences, presence-absence, or
  pseudo-absence, records for a given species

- pr_ab:

  character. Column name of "data" with presences, presence-absence, or
  pseudo-absence. Presences must be represented by 1 and absences by 0

- method:

  character. Vector with data partitioning method to be used. Usage
  part=c(method= 'kfold', folds='5'). Methods include:

  - kfold: Random partitioning into k-folds for cross-validation.
    'folds' refers to the number of folds for data partitioning, it
    assumes value \>=1. Usage method = c(method = "kfold", folds = 10).

  - rep_kfold: Random partitioning into repeated k-folds for
    cross-validation. Usage method = c(method = "rep_kfold", folds = 10,
    replicates=10). 'folds' refers to the number of folds for data
    partitioning, it assumes value \>=1. 'replicate' refers to the
    number of replicates, it assumes a value \>=1.

  - loocv: Leave-one-out cross-validation (a.k.a. Jackknife). It is a
    special case of k-fold cross validation where the number of
    partitions is equal to the number of records. Usage method =
    c(method = "loocv").

  - boot: Random bootstrap partitioning. Usage method=c(method='boot',
    replicates='2', proportion='0.7'). 'replicate' refers to the number
    of replicates, it assumes a value \>=1. 'proportion' refers to the
    proportion of occurrences used for model fitting, it assumes a value
    \>0 and \<=1. In this example proportion='0.7' mean that 70% of data
    will be used for model training, while 30% will be used for model
    testing.

## Value

A tibble object with information used in the 'data' argument and
additional columns named .part containing the partition groups. The
rep_kfold and boot method will return as many ".part" columns as
replicated defined. For the rest of the methods, a single .part column
is returned. For kfold, rep_kfold, and loocv partition methods, groups
are defined by integers. In contrast, for boot method, the partition
groups are defined by the characters 'train' and 'test'.

## References

- Fielding, A. H., & Bell, J. F. (1997). A review of methods for the
  assessment of prediction errors in conservation presence/absence
  models. Environmental Conservation, 24(1), 38-49.
  https://doi.org/10.1017/S0376892997000088

## See also

[`part_sblock`](https://sjevelazco.github.io/flexsdm/reference/part_sblock.md),
[`part_senv`](https://sjevelazco.github.io/flexsdm/reference/part_senv.md),
[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md),
[`sample_background`](https://sjevelazco.github.io/flexsdm/reference/sample_background.md)

## Examples

``` r
# \donttest{
require(dplyr)
data("abies")
abies$partition <- NULL
abies <- tibble(abies)

# K-fold method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 10)
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

# Repeated K-fold method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 10, replicates = 10)
)
abies2
#> # A tibble: 1,400 × 23
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
#> # ℹ 12 more variables: depth <dbl>, landform <fct>, .part1 <int>, .part2 <int>,
#> #   .part3 <int>, .part4 <int>, .part5 <int>, .part6 <int>, .part7 <int>,
#> #   .part8 <int>, .part9 <int>, .part10 <int>

# Leave-one-out cross-validation (loocv) method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "loocv")
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

# Bootstrap method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 50, proportion = 0.7)
)
abies2
#> # A tibble: 1,400 × 63
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
#> # ℹ 52 more variables: depth <dbl>, landform <fct>, .part1 <chr>, .part2 <chr>,
#> #   .part3 <chr>, .part4 <chr>, .part5 <chr>, .part6 <chr>, .part7 <chr>,
#> #   .part8 <chr>, .part9 <chr>, .part10 <chr>, .part11 <chr>, .part12 <chr>,
#> #   .part13 <chr>, .part14 <chr>, .part15 <chr>, .part16 <chr>, .part17 <chr>,
#> #   .part18 <chr>, .part19 <chr>, .part20 <chr>, .part21 <chr>, .part22 <chr>,
#> #   .part23 <chr>, .part24 <chr>, .part25 <chr>, .part26 <chr>, …
abies2$.part1 %>% table() # Note that for this method .partX columns have train and test words.
#> .
#>       test      train train-test 
#>        114        672        306 
# }
```
