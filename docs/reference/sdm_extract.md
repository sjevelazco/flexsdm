# Extract environmental data values from a spatial raster based on x and y coordinates

Extract environmental data values from a spatial raster based on x and y
coordinates

## Usage

``` r
sdm_extract(data, x, y, env_layer, variables = NULL, filter_na = TRUE)
```

## Arguments

- data:

  data.frame. Database with species presence, presence-absence, or
  pseudo-absence records with x and y coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- env_layer:

  SpatRaster. Raster or raster stack with environmental variables.

- variables:

  character. Vector with the variable names of predictor (environmental)
  variables Usage variables. = c("aet", "cwd", "tmin"). If no variable
  is specified, function will return data for all layers. Default NULL

- filter_na:

  logical. If filter_na = TRUE (default), the rows with NA values for
  any of the environmental variables are removed from the returned
  tibble.

## Value

A tibble that returns the original data base with additional columns for
the extracted environmental variables at each xy location from the
SpatRaster object used in 'env_layer'

## Examples

``` r
# \donttest{
require(terra)

# Load datasets
data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)

# Extract environmental data from somevar for all locations in spp
ex_spp <-
  sdm_extract(
    data = spp,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = NULL,
    filter_na = FALSE
  )

# Extract environmental for two variables and remove rows with NAs
ex_spp2 <-
  sdm_extract(
    data = spp,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = c("CFP_3", "CFP_4"),
    filter_na = TRUE
  )
#> 65 rows were excluded from database because NAs were found

ex_spp
#> # A tibble: 1,150 × 8
#>    species        x        y pr_ab CFP_1 CFP_2 CFP_3  CFP_4
#>    <chr>      <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
#>  1 sp1       -5541. -145138.     0 1160.  9.72  43.5  0.976
#>  2 sp1      -51981.   16322.     0  865.  8.27 138.   3.45 
#>  3 sp1     -269871.   69512.     1  746.  7.80 316.   4.67 
#>  4 sp1      -96261.  -32008.     0  980.  9.67  60.0  0.969
#>  5 sp1      269589. -566338.     0 1091. 12.0   54.6  1.18 
#>  6 sp1       29829. -328468.     0 1189.  7.84  40.2  1.09 
#>  7 sp1     -152691.  393782.     0  409.  2.06 176.  15.4  
#>  8 sp1     -195081.  253652.     0  830.  9.64 144.   8.27 
#>  9 sp1        -951. -277978.     0 1136. 10.2   47.4  0.875
#> 10 sp1      145929. -271498.     0 1031.  7.43  68.1  4.22 
#> # ℹ 1,140 more rows
ex_spp2
#> # A tibble: 1,085 × 6
#>    species        x        y pr_ab CFP_3  CFP_4
#>    <chr>      <dbl>    <dbl> <dbl> <dbl>  <dbl>
#>  1 sp1       -5541. -145138.     0  43.5  0.976
#>  2 sp1      -51981.   16322.     0 138.   3.45 
#>  3 sp1     -269871.   69512.     1 316.   4.67 
#>  4 sp1      -96261.  -32008.     0  60.0  0.969
#>  5 sp1      269589. -566338.     0  54.6  1.18 
#>  6 sp1       29829. -328468.     0  40.2  1.09 
#>  7 sp1     -152691.  393782.     0 176.  15.4  
#>  8 sp1     -195081.  253652.     0 144.   8.27 
#>  9 sp1        -951. -277978.     0  47.4  0.875
#> 10 sp1      145929. -271498.     0  68.1  4.22 
#> # ℹ 1,075 more rows
# }
```
