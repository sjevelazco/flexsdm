# A data set containing localities of Hesperocyparis stephensonii species in California, USA

A data set containing localities of Hesperocyparis stephensonii species
in California, USA

## Usage

``` r
hespero
```

## Format

A tibble object with 14 rows and 4 variables:

- ID:

  presences records ID

- x y:

  columns with coordinates in Albers Equal Area Conic coordinate system

- pr_ab:

  presence denoted by 1

## Examples

``` r
# \donttest{
require(dplyr)
data("hespero")
hespero
#> # A tibble: 21 × 4
#>       id       x        y pr_ab
#>    <int>   <dbl>    <dbl> <dbl>
#>  1     1 316923. -557843.     1
#>  2     2 317155. -559234.     1
#>  3     3 316960. -558186.     1
#>  4     4 314347. -559648.     1
#>  5     5 317348. -557349.     1
#>  6     6 316753. -559679.     1
#>  7     7 316777. -558644.     1
#>  8     8 317050. -559043.     1
#>  9     9 316655. -559928.     1
#> 10    10 316418. -567439.     1
#> # ℹ 11 more rows
# }
```
