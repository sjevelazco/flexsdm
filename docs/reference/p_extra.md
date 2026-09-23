# Graphical exploration of extrapolation or suitability pattern in the environmental and geographical space

Graphical exploration of extrapolation or suitability pattern in the
environmental and geographical space

## Usage

``` r
p_extra(
  training_data,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data,
  projection_data,
  predictors = NULL,
  geo_space = TRUE,
  geo_position = "right",
  prop_points = 0.2,
  maxcells = 1e+05,
  alpha_p = 0.5,
  color_p = "black",
  alpha_gradient = 0.5,
  color_gradient = c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E",
    "#30678D", "#3D4988", "#462777", "#440154"),
  theme = ggplot2::theme_classic()
)
```

## Arguments

- training_data:

  data.frame. Database with response (0,1) and predictor values used to
  fit a model.

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- pr_ab:

  character. Column name with species absence-presence,
  pseudo-absence-presence, or background-presence data (0,1).

- extra_suit_data:

  SpatRaster. Raster layer with extrapolation or suitability values.
  extra_suit_data must have the same resolution and extent than
  projection_data

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. projection_data must have the same resolution and extent
  than extra_suit_data

- predictors:

  character. Vector of predictor name(s) to calculate partial dependence
  plots. If NULL all predictors will be used. Default NULL.

- geo_space:

  logical. If TRUE will be produced a map. Default TRUE

- geo_position:

  character. Map position regarding plot of environmental space, right,
  left, bottom, or upper. Default "right"

- prop_points:

  numeric. Proportion of cells from extra_suit_data and projection_data
  to select for plotting. default. 0.5.

- maxcells:

  integer. Maximum number of cells used to plot in the geographical
  space. Default 100000

- alpha_p:

  numeric. a value between 0 to 1 to control transparency of
  presence-absence points. Lower values corresponding to more
  transparent colors. Default 0.5

- color_p:

  character. A vector with a color used to color presence-absence
  points. Default "black"

- alpha_gradient:

  numeric. a value between 0 to 1 to control transparency of projection
  data Lower values corresponding to more transparent colors. Default
  0.5

- color_gradient:

  character. A vector with colors used to color projection data. Default
  c( "#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E",
  "#30678D", "#3D4988", "#462777", "#440154")

- theme:

  ggplot2 theme. Default ggplot2::theme_classic()

## Value

a plot

## Examples

``` r
# \donttest{

require(dplyr)
require(terra)
require(ggplot2)

data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)
names(somevar) <- c("aet", "cwd", "tmx", "tmn")

spp$species %>% unique()
#> [1] "sp1" "sp2" "sp3"
sp <- spp %>%
  dplyr::filter(species == "sp2", pr_ab == 1) %>%
  dplyr::select(x, y, pr_ab)

# Calibration area based on some criterion such as dispersal ability
ca <- calib_area(sp,
  x = "x", y = "y",
  method = c("buffer", width = 50000), crs = crs(somevar)
)

plot(somevar[[1]])
points(sp)
plot(ca, add = TRUE)



# Sampling pseudo-absences
set.seed(10)
psa <- sample_pseudoabs(
  data = sp,
  x = "x",
  y = "y",
  n = nrow(sp) * 2,
  method = "random",
  rlayer = somevar,
  calibarea = ca
)

# Merge presences and abasences databases to get a complete calibration data
sp_pa <- dplyr::bind_rows(sp, psa)
sp_pa
#> # A tibble: 75 × 3
#>           x        y pr_ab
#>       <dbl>    <dbl> <dbl>
#>  1 -169701.   40352.     1
#>  2 -175101.  182642.     1
#>  3 -152151.  110552.     1
#>  4   34149. -168088.     1
#>  5  -88701. -248008.     1
#>  6  -63051.  -45508.     1
#>  7 -177801.  138632.     1
#>  8 -196431.  187772.     1
#>  9   54399. -125158.     1
#> 10  -35511.  -47128.     1
#> # ℹ 65 more rows

# Get environmental condition of calibration area
sp_pa_2 <- sdm_extract(data = sp_pa, x = "x", y = "y", env_layer = somevar)
sp_pa_2
#> # A tibble: 75 × 7
#>           x        y pr_ab   aet   cwd   tmx   tmn
#>       <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 -169701.   40352.     1  813.  9.02 128.  1.16 
#>  2 -175101.  182642.     1  745.  9.33 101.  4.52 
#>  3 -152151.  110552.     1  878.  9.64  93.8 2.22 
#>  4   34149. -168088.     1 1167.  9.90  45.5 1.20 
#>  5  -88701. -248008.     1  972.  5.15  90.8 0.698
#>  6  -63051.  -45508.     1 1023.  9.69  62.6 1.53 
#>  7 -177801.  138632.     1  867.  9.16  85.8 2.58 
#>  8 -196431.  187772.     1  895.  9.46 101.  4.21 
#>  9   54399. -125158.     1  993.  9.13 104.  2.66 
#> 10  -35511.  -47128.     1 1043.  9.44  66.7 1.51 
#> # ℹ 65 more rows

# Measure extrapolation based on calibration data (presence and pseudo-absences)
# using SHAPE metric
extr <-
  extra_eval(
    training_data = sp_pa_2,
    pr_ab = "pr_ab",
    projection_data = somevar,
    metric = "mahalanobis",
    univar_comb = FALSE,
    aggreg_factor = 1
  )
plot(extr)


## %######################################################%##
####            Explore extrapolation in the            ####
####        environmental and geographical space        ####
## %######################################################%##

p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05
)
#> Number of cell used to plot 3642 (5%)


p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr,
  projection_data = somevar,
  predictors = c("tmn", "cwd"),
  geo_space = TRUE,
  prop_points = 0.05
)
#> Number of cell used to plot 3642 (5%)


p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr,
  projection_data = somevar,
  predictors = c("cwd", "tmx", "aet"),
  geo_space = TRUE,
  geo_position = "left",
  prop_points = 0.05,
  color_p = "white",
  alpha_p = 0.5,
  alpha_gradient = 0.2,
  color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
  theme = ggplot2::theme_dark()
)
#> Number of cell used to plot 3642 (5%)


p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_p = "white",
  alpha_p = 0.5,
  alpha_gradient = 0.2,
  color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
  theme = ggplot2::theme_dark()
)
#> Number of cell used to plot 3642 (5%)


# Explore extrapolation only in the environmental space
p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr,
  projection_data = somevar,
  geo_space = FALSE,
  prop_points = 0.05,
  color_p = "black",
  color_gradient = c("#085CF8", "#65AF1E", "#F3CC1D", "#FC6A9B", "#D70500"),
  theme = ggplot2::theme_minimal()
)
#> Number of cell used to plot 3642 (5%)



## %######################################################%##
####                 Explore univariate                 ####
####          and combinatorial extrapolation           ####
## %######################################################%##
extr <-
  extra_eval(
    training_data = sp_pa_2,
    pr_ab = "pr_ab",
    projection_data = somevar,
    metric = "mahalanobis",
    univar_comb = TRUE,
    aggreg_factor = 1
  )

plot(extr)



p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = extr$uni_comb, # use  uni_comb layer
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_gradient = c("#B3DC2B", "#25818E")
)
#> Number of cell used to plot 3642 (5%)


## %######################################################%##
####           With p_extra also is possible            ####
####       to explore the patterns of suitability       ####
## %######################################################%##

sp_pa_2 <- part_random(
  data = sp_pa_2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

rf_m1 <- fit_raf(
  data = sp_pa_2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  partition = ".part",
  thr = c("max_sorensen")
)
#> Formula used for model fitting:
#> pr_ab ~ aet + cwd + tmx + tmn
#> Replica number: 1/1
#> Partition number: 1/5
#> Partition number: 2/5
#> Partition number: 3/5
#> Partition number: 4/5
#> Partition number: 5/5

suit <- sdm_predict(models = rf_m1, pred = somevar)
#> Predicting individual models
plot(suit$raf)

suit <- suit$raf

# Pasterns of suitability in geographical and environmental space
p_extra(
  training_data = sp_pa_2,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = suit,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
)
#> Number of cell used to plot 3642 (5%)


# Pasterns of suitability plotting as points only presences
p_extra(
  training_data = sp_pa_2 %>%
    dplyr::filter(pr_ab == 1),
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = suit,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
)
#> Number of cell used to plot 3642 (5%)


# Pasterns of suitability in the environmental space only
# and plotting as points only presences
p_extra(
  training_data = sp_pa_2 %>%
    dplyr::filter(pr_ab == 1),
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  extra_suit_data = suit,
  projection_data = somevar,
  geo_space = FALSE,
  prop_points = 0.05,
)
#> Number of cell used to plot 3642 (5%)

# }
```
