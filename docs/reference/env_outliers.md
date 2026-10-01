# Integration of outliers detection methods in environmental space

This function performs different methods for detecting outliers in
species distribution data based on the environmental conditions of
occurrences. Some methods need presence and absence data (e.g. Two-class
Support Vector Machine and Random Forest) while other only use presences
(e.g. Reverse Jackknife, Box-plot, and Random Forest outliers) . Outlier
detection can be a useful procedure in occurrence data cleaning (Chapman
2005, Liu et al., 2018).

## Usage

``` r
env_outliers(data, x, y, pr_ab, id, env_layer)
```

## Arguments

- data:

  data.frame or tibble with presence (or presence-absence) records, and
  coordinates

- x:

  character. Column name with longitude data.

- y:

  character. Column name with latitude data.

- pr_ab:

  character. Column name with presence and absence data (i.e. 1 and 0)

- id:

  character. Column name with row id. Each row (record) must have its
  own unique code.

- env_layer:

  SpatRaster. Raster with environmental variables

## Value

A tibble object with the same database used in 'data' argument and with
seven additional columns, where 1 and 0 denote that a presence was
detected or not as outliers

- .out_bxpt: outliers detected with Box-plot method

- .out_jack: outliers detected with Reverse Jackknife method

- .out_svm: outliers detected with Support Vector Machine method

- .out_rf: outliers detected with Random Forest method

- .out_rfout: outliers detected with Random Forest Outliers method

- .out_sum: frequency of a presences records was detected as outliers
  based on the previews methods (values between 0 and 6).

## Details

This function will apply outliers detection methods to occurrence data.
Box-plot and Reverse Jackknife method will test outliers for each
variable individually, if an occurrence behaves as an outlier for at
least one variable it will be highlighted as an outlier. If the user
uses only presence data, Support Vector Machine and Random Forest
Methods will not be performed. Support Vector Machine and Random Forest
are performed with default hyper-parameter values. In the case of a
species with \< 7 occurrences, the function will not perform any methods
(i.e. the additional columns will have 0 values); nonetheless, it will
return a tibble with the additional columns with 0 and 1. For further
information about these methods, see Chapman (2005), Liu et al. (2018),
and Velazco et al. (2022).

## References

- Chapman, A. D. (2005). Principles and methods of data cleaning:
  Primary Species and Species- Occurrence Data. version 1.0. Report for
  the Global Biodiversity Information Facility, Copenhagen. p72.
  http://www.gbif.org/document/80528

- Liu, C., White, M., & Newell, G. (2018). Detecting outliers in species
  distribution data. Journal of Biogeography, 45(1), 164 - 176.
  https://doi.org/10.1111/jbi.13122

- Velazco, S.J.E.; Bedrij, N.A.; Keller, H.A.; Rojas, J.L.; Ribeiro,
  B.R.; De Marco, P. (2022) Quantifying the role of protected areas for
  safeguarding the uses of biodiversity. Biological Conservation, xx(xx)
  xx-xx. https://doi.org/10.1016/j.biocon.2022.109525

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)
require(ggplot2)
#> Loading required package: ggplot2
#> Warning: package 'ggplot2' was built under R version 4.5.3

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Species occurrences
data("spp")
spp
#> # A tibble: 1,150 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1       -5541. -145138.     0
#>  2 sp1      -51981.   16322.     0
#>  3 sp1     -269871.   69512.     1
#>  4 sp1      -96261.  -32008.     0
#>  5 sp1      269589. -566338.     0
#>  6 sp1       29829. -328468.     0
#>  7 sp1     -152691.  393782.     0
#>  8 sp1     -195081.  253652.     0
#>  9 sp1        -951. -277978.     0
#> 10 sp1      145929. -271498.     0
#> # ℹ 1,140 more rows
spp1 <- spp %>% dplyr::filter(species == "sp1")

somevar[[1]] %>% plot()
points(spp1 %>% filter(pr_ab == 1) %>% select(x, y), col = "blue", pch = 19)
points(spp1 %>% filter(pr_ab == 0) %>% select(x, y), col = "red", cex = 0.5)


spp1 <- spp1 %>% mutate(idd = 1:nrow(spp1))

# Detect outliers
outs_1 <- env_outliers(
  data = spp1,
  pr_ab = "pr_ab",
  x = "x",
  y = "y",
  id = "idd",
  env_layer = somevar
)
#> 55 rows were excluded from database because NAs were found

# How many outliers were detected by different methods?
out_pa <- outs_1 %>%
  dplyr::select(starts_with("."), -.out_sum) %>%
  apply(., 2, function(x) sum(x, na.rm = TRUE))
out_pa
#>  .out_bxpt  .out_jack   .out_svm    .out_rf .out_rfout 
#>         19          0         12         12         12 

# How many outliers were detected by the sum of different methods?
outs_1 %>%
  dplyr::group_by(.out_sum) %>%
  dplyr::count()
#> # A tibble: 5 × 2
#> # Groups:   .out_sum [5]
#>   .out_sum     n
#>      <dbl> <int>
#> 1        0   903
#> 2        1    31
#> 3        2     9
#> 4        3     2
#> 5       NA    55

# Let explor where are locate records highlighted as outliers
outs_1 %>%
  dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(col = factor(.out_sum))) +
  facet_wrap(. ~ factor(.out_sum))


# Detect outliers only with presences
outs_2 <- env_outliers(
  data = spp1 %>% dplyr::filter(pr_ab == 1),
  pr_ab = "pr_ab",
  x = "x",
  y = "y",
  id = "idd",
  env_layer = somevar
)
#> 12 rows were excluded from database because NAs were found

# How many outliers were detected by different methods
out_p <- outs_2 %>%
  dplyr::select(starts_with("."), -.out_sum) %>%
  apply(., 2, function(x) sum(x, na.rm = TRUE))

# How many outliers were detected by the sum of different methods?
outs_2 %>%
  dplyr::group_by(.out_sum) %>%
  dplyr::count()
#> # A tibble: 4 × 2
#> # Groups:   .out_sum [4]
#>   .out_sum     n
#>      <dbl> <int>
#> 1        0   208
#> 2        1    29
#> 3        2     1
#> 4       NA    12

# Let explor where are locate records highlighted as outliers
outs_2 %>%
  dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(col = factor(.out_sum))) +
  facet_wrap(. ~ factor(.out_sum))



# Comparison of function outputs when using it with
# presences-absences or only presences data.

bind_rows(out_p, out_pa)
#> # A tibble: 2 × 5
#>   .out_bxpt .out_jack .out_svm .out_rf .out_rfout
#>       <dbl>     <dbl>    <dbl>   <dbl>      <dbl>
#> 1        19         0        0       0         12
#> 2        19         0       12      12         12
# Because the second case only were used presences, outliers methods
# based in Random Forest (.out_rf) and Support Vector Machines (.out_svm)
# were not performed.
# }
```
