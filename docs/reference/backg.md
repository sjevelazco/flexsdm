# A data set containing environmental conditions of background points

A data set containing environmental conditions of background points

## Usage

``` r
backg
```

## Format

A tibble object with 5000 rows and 10 variables:

- pr_ab:

  background point denoted by 0

- x y:

  columns with geographical coordinates

- from column aet to landform:

  columns with values of environmental variables at coordinate locations

## Examples

``` r
# \donttest{
require(dplyr)
data("backg")
backg
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
# }
```
