# A data set containing localities and environmental condition of an Abies (fir tree) species in California, USA

A data set containing localities and environmental condition of an Abies
(fir tree) species in California, USA

## Usage

``` r
abies
```

## Format

A tibble object with 5000 rows and 10 variables:

- ID:

  presences and absences records ID

- pr_ab:

  presence and absences denoted by 1 and 0 respectively

- x y:

  columns with coordinates in Albers Equal Area Conic coordinate system

- from column aet to landform:

  columns with values for environmental variables at each locality

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
data("abies")
abies
} # }
```
