skip_on_cran() # slow tests; they run when NOT_CRAN=true

require(terra)
require(dplyr)

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Species occurrences
data("spp")
spp
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

## %######################################################%##
####                  Cellsize method                   ####
## %######################################################%##
# Using cellsize method with several factors. Thinning is random and occfilt_geo()
# does not set a seed, so each factor is run with the same seed to get reproducible results
cellsize_factors <- c(1, 4, 8, 12, 16, 20)
filtered_occ <- lapply(cellsize_factors, function(f) {
  set.seed(1)
  occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("cellsize", factor = f),
    prj = crs(somevar)
  )
})
names(filtered_occ) <- cellsize_factors


# Select filtered occurrences based on
# number of records and spatial autocorrelation

test_that("test occfilt_select filter_prop = FALSE", {
  occ_selected <- occfilt_select(
    occ_list = filtered_occ,
    x = "x",
    y = "y",
    env_layer = somevar,
    filter_prop = FALSE
  )

  expect_equal(class(occ_selected)[1], "tbl_df")
  expect_equal(nrow(occ_selected), 118)
})


test_that("test occfilt_select with filter_prop = TRUE", {
  skip_on_cran()
  occ_selected <- occfilt_select(
    occ_list = filtered_occ,
    x = "x",
    y = "y",
    env_layer = somevar,
    filter_prop = TRUE
  )

  expect_equal(class(occ_selected), "list")
  expect_equal(length(occ_selected), 2)
  expect_equal(nrow(occ_selected$occ), 118)
  expect_equal(nrow(occ_selected$filter_prop), 6)
})

test_that("expect error", {
  skip_on_cran()
  expect_error(occfilt_select(
    occ_list = filtered_occ[[1]],
    x = "x",
    y = "y",
    env_layer = somevar,
    filter_prop = TRUE
  ))
})
