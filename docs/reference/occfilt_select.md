# Select filtered occurrences

Select filtered occurrences based on number of records and spatial
autocorrelation (see details)

## Usage

``` r
occfilt_select(occ_list, x, y, env_layer, filter_prop = FALSE)
```

## Arguments

- occ_list:

  list. A list with filtered specie occurrences testing several values
  (see
  [`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md)
  and
  [`occfilt_geo`](https://sjevelazco.github.io/flexsdm/reference/occfilt_geo.md))

- x:

  character. Column name with longitude data

- y:

  character. Column name with latitude data

- env_layer:

  SpatRaster. Raster variables that will be used to fit the model.
  Factor variables will be removed.

- filter_prop:

  logical. If TRUE, the function will return a list with the filtered
  occurrences and a tibble with the spatial autocorrelation and number
  of occurrence values

## Value

If filter_prop = FALSE, a tibble with selected filtered occurrences. If
filter_prop = TRUE, a list with following objects:

- A tibble with selected filtered occurrences

- A tibble with filter properties with columns:

  - filt_value: values used for filtering, the value with an asterisk
    will denote the one selected

  - n_records: number of occurrence

  - mean_autocorr: mean spatial autocorrelation.

  - the remaining columns have the spatial autocorrelation values for
    each variable.

## Details

The function implement the approach used in Velazco et al. (2020) which
consists in calculating for each filtered dataset:

- 1- the number of occurrence.

- 2- the spatial autocorrelation based on Morans'I for each variable

- 3- the mean spatial autocorrelation among variables

Then function will select those dataset with average spatial
autocorrelation lower than the mean of all dataset, and from this subset
will select the one with the highest number occurrences.

If use occfilt_select cite Velazco et al. (2020) as reference.

## References

- Velazco, S. J. E., Svenning, J-C., Ribeiro, B. R., & Laureto, L. M. O.
  (2020). On opportunities and threats to conserve the phylogenetic
  diversity of Neotropical palms. Diversity and Distributions, 27,
  512–523. https://doi.org/10.1111/ddi.13215

## See also

[`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md),
[`occfilt_geo`](https://sjevelazco.github.io/flexsdm/reference/occfilt_geo.md)

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

plot(somevar)


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
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

## %######################################################%##
####                  Cellsize method                   ####
## %######################################################%##
# Using cellsize method
filtered_occ <- occfilt_geo(
  data = spp1,
  x = "x",
  y = "y",
  env_layer = somevar,
  method = c("cellsize", factor = c(1, 4, 8, 12, 16, 20)),
  prj = crs(somevar)
)
#> Extracting values from raster ... 
#> 16 records were removed because they have NAs for some variables
#> Number of unfiltered records: 234
#> Factor: x1
#> Distance threshold (km): 1.539
#> Number of filtered records: 233
#> Factor: x4
#> Distance threshold (km): 6.157
#> Number of filtered records: 203
#> Factor: x8
#> Distance threshold (km): 12.313
#> Number of filtered records: 156
#> Factor: x12
#> Distance threshold (km): 18.47
#> Number of filtered records: 118
#> Factor: x16
#> Distance threshold (km): 24.626
#> Number of filtered records: 96
#> Factor: x20
#> Distance threshold (km): 30.783
#> Number of filtered records: 78

filtered_occ
#> $`1`
#> # A tibble: 233 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -126231.  196142.     1
#>  4 sp1       91659. -156748.     1
#>  5 sp1     -210471.  326282.     1
#>  6 sp1     -140541.  284972.     1
#>  7 sp1     -217491.   65732.     1
#>  8 sp1     -178611.  225032.     1
#>  9 sp1      -92481.  155642.     1
#> 10 sp1     -367611.  266072.     1
#> # ℹ 223 more rows
#> 
#> $`4`
#> # A tibble: 203 × 4
#>    species        x       y pr_ab
#>    <chr>      <dbl>   <dbl> <dbl>
#>  1 sp1     -269871.  69512.     1
#>  2 sp1     -149991. 267962.     1
#>  3 sp1     -126231. 196142.     1
#>  4 sp1     -210471. 326282.     1
#>  5 sp1     -140541. 284972.     1
#>  6 sp1     -217491.  65732.     1
#>  7 sp1     -178611. 225032.     1
#>  8 sp1      -92481. 155642.     1
#>  9 sp1     -367611. 266072.     1
#> 10 sp1     -109491.  46292.     1
#> # ℹ 193 more rows
#> 
#> $`8`
#> # A tibble: 156 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1       91659. -156748.     1
#>  4 sp1     -210471.  326282.     1
#>  5 sp1     -140541.  284972.     1
#>  6 sp1     -217491.   65732.     1
#>  7 sp1      -92481.  155642.     1
#>  8 sp1     -367611.  266072.     1
#>  9 sp1     -184551.  425372.     1
#> 10 sp1     -260151.   57632.     1
#> # ℹ 146 more rows
#> 
#> $`12`
#> # A tibble: 118 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -140541.  284972.     1
#>  4 sp1     -178611.  225032.     1
#>  5 sp1      -34431.  212072.     1
#>  6 sp1     -331701.  304412.     1
#>  7 sp1     -215871.  -51178.     1
#>  8 sp1      -60351. -307138.     1
#>  9 sp1     -331161.  387842.     1
#> 10 sp1      -88971.  137282.     1
#> # ℹ 108 more rows
#> 
#> $`16`
#> # A tibble: 96 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -210471.  326282.     1
#>  2 sp1     -313611.  400802.     1
#>  3 sp1     -202911.  -64138.     1
#>  4 sp1      -34431.  212072.     1
#>  5 sp1     -172941. -113008.     1
#>  6 sp1      -60351. -307138.     1
#>  7 sp1     -304971.  219632.     1
#>  8 sp1     -195351.   85982.     1
#>  9 sp1      -23091.   56012.     1
#> 10 sp1     -114621.   54122.     1
#> # ℹ 86 more rows
#> 
#> $`20`
#> # A tibble: 78 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -126231.  196142.     1
#>  2 sp1     -367611.  266072.     1
#>  3 sp1     -109491.   46292.     1
#>  4 sp1      -34431.  212072.     1
#>  5 sp1     -143781.  233942.     1
#>  6 sp1     -215871.  -51178.     1
#>  7 sp1      -60351. -307138.     1
#>  8 sp1     -306321.  134852.     1
#>  9 sp1     -195351.   85982.     1
#> 10 sp1     -328731.  405392.     1
#> # ℹ 68 more rows
#> 

# Select filtered occurrences based on
# number of records and spatial autocorrelation
occ_selected <- occfilt_select(
  occ_list = filtered_occ,
  x = "x",
  y = "y",
  env_layer = somevar,
  filter_prop = FALSE
)
#> Dataset with filtered value 12 was selected
occ_selected
#> # A tibble: 118 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -140541.  284972.     1
#>  4 sp1     -178611.  225032.     1
#>  5 sp1      -34431.  212072.     1
#>  6 sp1     -331701.  304412.     1
#>  7 sp1     -215871.  -51178.     1
#>  8 sp1      -60351. -307138.     1
#>  9 sp1     -331161.  387842.     1
#> 10 sp1      -88971.  137282.     1
#> # ℹ 108 more rows

occ_selected <- occfilt_select(
  occ_list = filtered_occ,
  x = "x",
  y = "y",
  env_layer = somevar,
  filter_prop = TRUE
)
#> Dataset with filtered value 12 was selected
occ_selected$occ
#> # A tibble: 118 × 4
#>    species        x        y pr_ab
#>    <chr>      <dbl>    <dbl> <dbl>
#>  1 sp1     -269871.   69512.     1
#>  2 sp1     -149991.  267962.     1
#>  3 sp1     -140541.  284972.     1
#>  4 sp1     -178611.  225032.     1
#>  5 sp1      -34431.  212072.     1
#>  6 sp1     -331701.  304412.     1
#>  7 sp1     -215871.  -51178.     1
#>  8 sp1      -60351. -307138.     1
#>  9 sp1     -331161.  387842.     1
#> 10 sp1      -88971.  137282.     1
#> # ℹ 108 more rows

occ_selected$filter_prop
#>   filt_value mean_autocorr n_records     CFP_1      CFP_2     CFP_3     CFP_4
#> 1          1     0.3219265       233 0.3534014 0.20178498 0.3352363 0.3972832
#> 2          4     0.2860045       203 0.3181670 0.15331928 0.3015718 0.3709598
#> 3          8     0.2602486       156 0.3055744 0.12870505 0.2635844 0.3431307
#> 4       * 12     0.2291901       118 0.2568785 0.10808796 0.2311247 0.3206691
#> 5         16     0.2136998        96 0.2664798 0.07228080 0.2457754 0.2702633
#> 6         20     0.2045660        78 0.2539334 0.07758901 0.1968988 0.2898429
# }
```
