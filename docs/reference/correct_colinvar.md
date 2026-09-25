# Collinearity reduction of predictor variables

Collinearity reduction of predictor variables

## Usage

``` r
correct_colinvar(
  env_layer,
  method,
  proj = NULL,
  save_proj = NULL,
  restric_to_region = NULL,
  restric_pca_proj = FALSE,
  maxcell = NULL,
  based_on_points = FALSE,
  data = NULL,
  x = NULL,
  y = NULL
)
```

## Arguments

- env_layer:

  SpatRaster An object of class SpatRaster containing the predictors.
  This function does not allow categorical variables

- method:

  character. Collinearity reduction method. It is necessary to provide a
  vector for this argument. The next methods are implemented:

  - pearson: Highlights correlated variables according to Pearson
    correlation. A threshold of maximum correlation must be specified.
    Otherwise, a threshold of 0.7 is defined as default. Usage method =
    c('pearson', th='0.7').

  - vif: Select variables by Variance Inflation Factor, a threshold can
    be specified by user. Otherwise, a threshold of 10 is defined as
    default.Usage method = c('vif', th = '10').

  - pca: Perform a Principal Component Analysis and use the principal
    components as the new predictors. The selected components account
    for 95% of the whole variation in the system. Usage method =
    c('pca').

  - fa: Perform a Factorial Analysis and select, from the original
    predictors, the number of factors is defined by Broken-Stick and
    variables with the highest correlation to the factors are selected.
    Usage method = c('fa').

- proj:

  character. Only used for pca method. Path to a folder that contains
  sub-folders for the different projection scenarios. Variables names
  must have the same names as in the raster used in env_layer argument.
  Usage proj = "C:/User/Desktop/Projections" (see in Details more about
  the use of this argument)

- save_proj:

  character. Directory to save PCA projection. Default NULL.

- restric_to_region:

  SpatVector. Area used to restrict cells of env_layer at moment to
  perform collinearity reduction. Default: NULL.

- restric_pca_proj:

  logical. Area used to restrict geographically PCA projection within
  SpatVector used in restric_to_region. Only use for PCA analysis.
  Default: FALSE.

- maxcell:

  numeric. Number of raster cells to be randomly sampled. Taking a
  sample could be useful to reduce memory usage for large rasters. If
  NULL, the function will use all raster cells. Default NULL. Usage
  maxcell = 50000.

- based_on_points:

  logical. If TRUE, collinearity reduction method will be based on
  species points data (i.e., presences, and absences, pseudo-absences or
  background points). If TRUE, data, x and y arguments must be provided.
  Default FALSE.

- data:

  tibble or data.frame. Database with species data used to model (i.e.,
  presence + absence, or presence + pseudo-absence + background points)
  with x and y coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

## Value

\#' If 'pearson', returns a list with the following elements:

- cor_table: a matrix object with pairwise correlation values of the
  environmental variables

- cor_variables: a list object with the same length of the number of
  environmental values containing the pairwise relations that exceeded
  the correlation threshold for each one of the environmental variables

If 'vif' method, returns a list with the following elements:

- env_layer: a SpatRaster object with selected environmental variables

- removed_variables: a character vector with removed environmental
  variables

- vif_table: a data frame with VIF values for all environmental
  variables

If 'pca' method, returns a list with the following elements:

- env_layer: SpatRaster with scores of selected principal component (PC)
  that sum up 95% of the whole variation or original environmental
  variables

- coefficients: a matrix with the coefficient of principal component
  (PC) for predictors

- cumulative_variance: a tibble with the cumulative variance explained
  in selected principal component (PC)

If 'fa' method, returns a list with the following elements:

- env_layer: SpatRaster with scores of selected variables due to
  correlation to factors.

- number_factors: number of factors selected according to the
  Broken-Stick criteria,

- removed_variables: removed variables,

- uniqueness: uniqueness of each environmental variable according to the
  factorial analysis,

- loadings: environmental variables loadings in each of the chosen
  factors

## Details

In the case of having environmental variables for the current conditions
and other time periods (future or present), it is recommended to perform
the PCA analysis with the current environmental condition and project
the PCA for the other time periods. To do so, it is necessary to use
“proj” argument. Path to a folder (e.g., projections) that contains
sub-folders for the different projection scenarios (e.g., years and
emissions). Within each sub-folder must be stored single or multiband
rasters with the environmental variables.

For example:

C:/Users/my_pc/projections/\
├── MRIESM_2050_ssp126\
│ └── var1.tif\
│ └── var2.tif\
│ └── var3.tif\
├── MRIESM_2080_ssp585\
│ └── var1.tif\
│ └── var2.tif\
│ └── var3.tif\
├── UKESM_2050_ssp370\
│ └── var1.tif\
│ └── var2.tif\
│ └── var3.tif

If pca method is run with time projections, correct_colinvar function
will create the Projection_PCA (the exact path is in the path object
returned by the function) with the same system of sub-folders and
multiband raster with the principal components (pcs.tif)

C:/Users/my_pc/Projection_PCA/\
├── MRIESM_2050_ssp126\
│ └── pcs.tif \# a multiband tif with principal components\
├── MRIESM_2080_ssp585\
│ └── pcs.tif\
├── UKESM_2050_ssp370\
│ └── pcs.tif

Perform collinearity reduction based on points

Evaluating collinearity based on all environmental conditions of a
calibration area or study area could yield different results than
evaluating collinearity based on points used to construct the models. If
you want to perform collinearity reduction based on species points data,
it is strongly recommended to use all the point data used for modeling
(i.e., presence + absence or presence + pseudo-absence/background
points).

## Examples

``` r
# \donttest{
require(terra)
require(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Perform pearson collinearity control
var <- correct_colinvar(env_layer = somevar, method = c("pearson", th = "0.7"))
var$cor_table
#>           CFP_1     CFP_2     CFP_3     CFP_4
#> CFP_1 0.0000000 0.7689893 0.7924813 0.7845401
#> CFP_2 0.7689893 0.0000000 0.4168956 0.5881831
#> CFP_3 0.7924813 0.4168956 0.0000000 0.7323259
#> CFP_4 0.7845401 0.5881831 0.7323259 0.0000000
var$cor_variables
#> $CFP_1
#> [1] "CFP_2" "CFP_3" "CFP_4"
#> 
#> $CFP_2
#> [1] "CFP_1"
#> 
#> $CFP_3
#> [1] "CFP_1" "CFP_4"
#> 
#> $CFP_4
#> [1] "CFP_1" "CFP_3"
#> 

# For all correct_colinvar methods it is possible to take a sample or raster to reduce memory
var <- correct_colinvar(env_layer = somevar, method = c("pearson", th = "0.7"), maxcell = 10000)
var$cor_table
#>           CFP_1     CFP_2     CFP_3     CFP_4
#> CFP_1 0.0000000 0.7653187 0.7945362 0.7855346
#> CFP_2 0.7653187 0.0000000 0.4155557 0.5885450
#> CFP_3 0.7945362 0.4155557 0.0000000 0.7303327
#> CFP_4 0.7855346 0.5885450 0.7303327 0.0000000
var$cor_variables
#> $CFP_1
#> [1] "CFP_2" "CFP_3" "CFP_4"
#> 
#> $CFP_2
#> [1] "CFP_1"
#> 
#> $CFP_3
#> [1] "CFP_1" "CFP_4"
#> 
#> $CFP_4
#> [1] "CFP_1" "CFP_3"
#> 

# Perform vif collinearity control
var <- correct_colinvar(env_layer = somevar, method = c("vif", th = "8"))
var$env_layer
#> class       : SpatRaster
#> size        : 558, 394, 4  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source      : somevar.tif
#> names       :       CFP_1,     CFP_2,      CFP_3,     CFP_4
#> min values  :           0,  -9.39489,  22.446851,  0.259143
#> max values  : 1357.864694, 14.200473, 614.691254, 64.374759
var$removed_variables
#> NULL
var$vif_table
#> # A tibble: 4 × 2
#>   Variables   VIF
#>   <chr>     <dbl>
#> 1 CFP_1      7.62
#> 2 CFP_2      3.29
#> 3 CFP_3      3.95
#> 4 CFP_4      2.89

# Perform pca collinearity control
var <- correct_colinvar(env_layer = somevar, method = c("pca"))
plot(var$env_layer)

var$env_layer
#> class       : SpatRaster
#> size        : 558, 394, 3  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> names       :       PC1,       PC2,       PC3
#> min values  : -8.453273, -4.260147, -1.525085
#> max values  :  2.827164,  3.337545,  4.342864
var$coefficients
#> # A tibble: 4 × 5
#>   variable    PC1     PC2    PC3     PC4
#>   <chr>     <dbl>   <dbl>  <dbl>   <dbl>
#> 1 CFP_1     0.550 -0.0722  0.296 -0.778 
#> 2 CFP_2     0.450 -0.777   0.103  0.429 
#> 3 CFP_3    -0.485 -0.594  -0.450 -0.459 
#> 4 CFP_4    -0.511 -0.198   0.836 -0.0241
var$cumulative_variance
#> # A tibble: 4 × 2
#>      PC  cvar
#>   <int> <dbl>
#> 1     1 0.764
#> 2     2 0.915
#> 3     3 0.979
#> 4     4 1    


# Perform pca collinearity control with different projections
## Below will be created a set of folders to simulate the structure of the
## directory where environmental variables are stored for different scenarios
dir_sc <- file.path(tempdir(), "projections")
dir.create(dir_sc)
dir_sc <- file.path(dir_sc, c("scenario_1", "scenario_2"))
sapply(dir_sc, dir.create)
#> /var/folders/0q/hy7d1sjn455f1kdjtvcpjg7h0000gn/T//RtmpKeEw59/projections/scenario_1 
#>                                                                                TRUE 
#> /var/folders/0q/hy7d1sjn455f1kdjtvcpjg7h0000gn/T//RtmpKeEw59/projections/scenario_2 
#>                                                                                TRUE 

somevar <-
  system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

terra::writeRaster(somevar, file.path(dir_sc[1], "somevar.tif"), overwrite = TRUE)
terra::writeRaster(somevar, file.path(dir_sc[2], "somevar.tif"), overwrite = TRUE)

## Perform pca with projections
dir_w_proj <- dirname(dir_sc[1])
dir_w_proj
#> [1] "/var/folders/0q/hy7d1sjn455f1kdjtvcpjg7h0000gn/T//RtmpKeEw59/projections"
var <- correct_colinvar(env_layer = somevar, method = "pca", proj = dir_w_proj)
var$env_layer
#> class       : SpatRaster
#> size        : 558, 394, 3  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> names       :       PC1,       PC2,       PC3
#> min values  : -8.453273, -4.260147, -1.525085
#> max values  :  2.827164,  3.337545,  4.342864
var$coefficients
#> # A tibble: 4 × 5
#>   variable    PC1     PC2    PC3     PC4
#>   <chr>     <dbl>   <dbl>  <dbl>   <dbl>
#> 1 CFP_1     0.550 -0.0722  0.296 -0.778 
#> 2 CFP_2     0.450 -0.777   0.103  0.429 
#> 3 CFP_3    -0.485 -0.594  -0.450 -0.459 
#> 4 CFP_4    -0.511 -0.198   0.836 -0.0241
var$cumulative_variance
#> # A tibble: 4 × 2
#>      PC  cvar
#>   <int> <dbl>
#> 1     1 0.764
#> 2     2 0.915
#> 3     3 0.979
#> 4     4 1    
var$proj
#> [1] "/var/folders/0q/hy7d1sjn455f1kdjtvcpjg7h0000gn/T//RtmpKeEw59/Projection_PCA"


# Perform fa colinearity control
var <- correct_colinvar(env_layer = somevar, method = c("fa"))
var$env_layer
#> class       : SpatRaster
#> size        : 558, 394, 1  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source      : somevar.tif
#> name        :       CFP_1
#> min value   :           0
#> max value   : 1357.864694
var$number_factors
#> [1] 1
var$removed_variables
#> [1] "CFP_2" "CFP_3" "CFP_4"
var$uniqueness
#>     CFP_1     CFP_2     CFP_3     CFP_4 
#> 0.0010000 0.4087429 0.3715653 0.3835564 
var$loadings
#> # A tibble: 4 × 2
#>   Variable Factor_1
#>   <chr>       <dbl>
#> 1 CFP_1       1.000
#> 2 CFP_2       0.769
#> 3 CFP_3      -0.793
#> 4 CFP_4      -0.785

## %######################################################%##
#                                                          #
####            Other option to perform PCA             ####
####      considering cell restricted to a region       ####
#                                                          #
## %######################################################%##
data("abies")

# Define a calibration area
abies2 <- abies %>%
  dplyr::select(x, y, pr_ab) %>%
  dplyr::filter(pr_ab == 1)

plot(somevar[[1]])
points(abies2[-3])
ca <- calib_area(abies2, x = "x", y = "y", method = c("mcp"), crs = crs(somevar))
plot(ca, add = TRUE)


# Full geographical range to perform PCA
pca_fr <- correct_colinvar(
  env_layer = somevar,
  method = c("pca"),
  maxcell = NULL,
  restric_to_region = NULL,
  restric_pca_proj = FALSE
)

# Perform PCA only with cell delimited by polygon used in restric_to_region
pca_rr <- correct_colinvar(
  env_layer = somevar,
  method = c("pca"),
  maxcell = NULL,
  restric_to_region = ca,
  restric_pca_proj = FALSE
)

# Perform and predicted PCA only with cell delimited by polygon used in restric_to_region
pca_rrp <- correct_colinvar(
  env_layer = somevar,
  method = c("pca"),
  maxcell = NULL,
  restric_to_region = ca,
  restric_pca_proj = TRUE
)

plot(pca_fr$env_layer) # PCA with all cells

plot(pca_rr$env_layer) # PCA with calibration area cell but predicted for entire region

plot(pca_rrp$env_layer) # PCA performed and predicted for cells within calibration area (ca)



## %######################################################%##
#                                                          #
####       Use correct_colinvar with points data        ####
#                                                          #
## %######################################################%##
data("abies")

# Presence-absence database
abies2 <- abies %>%
  dplyr::select(x, y, pr_ab)

# Perform collinearity control
# Pearson
correct_colinvar(
  env_layer = somevar,
  method = c("pearson", th = "0.6"),
  based_on_points = TRUE,
  data = abies2,
  x = "x",
  y = "y"
)
#> 60 rows were excluded from database because NAs were found
#> $cor_table
#>           CFP_1      CFP_2     CFP_3      CFP_4
#> CFP_1 0.0000000 0.65011233 0.4708649 0.55696032
#> CFP_2 0.6501123 0.00000000 0.2528434 0.09606692
#> CFP_3 0.4708649 0.25284342 0.0000000 0.84185543
#> CFP_4 0.5569603 0.09606692 0.8418554 0.00000000
#> 
#> $cor_variables
#> $cor_variables$CFP_1
#> [1] "CFP_2"
#> 
#> $cor_variables$CFP_2
#> [1] "CFP_1"
#> 
#> $cor_variables$CFP_3
#> [1] "CFP_4"
#> 
#> $cor_variables$CFP_4
#> [1] "CFP_3"
#> 
#> 

# VIF
correct_colinvar(
  env_layer = somevar,
  method = c("vif", th = "8"),
  based_on_points = TRUE,
  data = abies2,
  x = "x",
  y = "y"
)
#> 60 rows were excluded from database because NAs were found
#> $env_layer
#> class       : SpatRaster
#> size        : 558, 394, 4  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source      : somevar.tif
#> names       :       CFP_1,     CFP_2,      CFP_3,     CFP_4
#> min values  :           0,  -9.39489,  22.446851,  0.259143
#> max values  : 1357.864694, 14.200473, 614.691254, 64.374759
#> 
#> $removed_variables
#> NULL
#> 
#> $vif_table
#> # A tibble: 4 × 2
#>   Variables   VIF
#>   <chr>     <dbl>
#> 1 CFP_1      7.62
#> 2 CFP_2      5.91
#> 3 CFP_3      5.54
#> 4 CFP_4      4.03
#> 

# PCA
correct_colinvar(
  env_layer = somevar,
  method = c("pca"),
  based_on_points = TRUE,
  data = abies2,
  x = "x",
  y = "y"
)
#> 60 rows were excluded from database because NAs were found
#> $env_layer
#> class       : SpatRaster
#> size        : 558, 394, 3  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> names       :       PC1,       PC2,       PC3
#> min values  : -5.059408, -3.126523, -2.134568
#> max values  :  6.396155,  3.776424,  1.788421
#> 
#> $coefficients
#> # A tibble: 4 × 5
#>   variable    PC1   PC2     PC3     PC4
#>   <chr>     <dbl> <dbl>   <dbl>   <dbl>
#> 1 CFP_1    -0.536 0.453 -0.237   0.672 
#> 2 CFP_2    -0.113 0.792  0.0902 -0.592 
#> 3 CFP_3     0.576 0.343  0.597   0.440 
#> 4 CFP_4     0.606 0.222 -0.761   0.0652
#> 
#> $cumulative_variance
#> # A tibble: 4 × 2
#>      PC  cvar
#>   <int> <dbl>
#> 1     1 0.569
#> 2     2 0.945
#> 3     3 0.984
#> 4     4 1    
#> 

# FA
correct_colinvar(
  env_layer = somevar,
  method = "fa",
  based_on_points = TRUE,
  data = abies2,
  x = "x",
  y = "y"
)
#> 60 rows were excluded from database because NAs were found
#> Factorial analysis could not be performed because 2 factors are too many for 4 variables
#> It will tested with 1 factors
#> $env_layer
#> class       : SpatRaster
#> size        : 558, 394, 1  (nrow, ncol, nlyr)
#> resolution  : 1890, 1890  (x, y)
#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
#> source      : somevar.tif
#> name        :     CFP_4
#> min value   :  0.259143
#> max value   : 64.374759
#> 
#> $number_factors
#> [1] 1
#> 
#> $removed_variables
#> [1] "CFP_1" "CFP_2" "CFP_3"
#> 
#> $uniqueness
#>     CFP_1     CFP_2     CFP_3     CFP_4 
#> 0.6895565 0.9907748 0.2905398 0.0010000 
#> 
#> $loadings
#> # A tibble: 4 × 2
#>   Variable Factor_1
#>   <chr>       <dbl>
#> 1 CFP_1     -0.557 
#> 2 CFP_2      0.0960
#> 3 CFP_3      0.842 
#> 4 CFP_4      0.999 
#> 
# }
```
