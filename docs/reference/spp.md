# A data set containing presences and absences of three virtual species

A data set containing presences and absences of three virtual species

## Usage

``` r
spp
```

## Format

A tibble with 1150 rows and 3 variables:

- species:

  virtual species names

- x:

  longitude of species occurrences

- y:

  latitude of species occurrences

- pr_ab:

  presences and absences denoted by 1 and 0 respectively

## Examples

``` r
# \donttest{
require(dplyr)
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
# }
```
