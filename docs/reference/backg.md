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
if (FALSE) { # \dontrun{
require(dplyr)
data("backg")
backg
} # }
```
