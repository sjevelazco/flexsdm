require(terra)
require(dplyr)

# Environmental variables
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)

# Species occurrences
data("spp")
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)


test_that("occfilt_geo 'moran' method", {
  # Using Moran method
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("moran"),
    prj = crs(somevar)
  )
  expect_equal(nrow(filtered), 4)
})

test_that("occfilt_geo 'cellsize' method for different values", {
  # Test for different values
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("moran", c(0.05, 0.05, 0.2)),
    prj = crs(somevar)
  )
  expect_equal(class(filtered), "list")
  expect_equal(length(filtered), 3)
})



test_that("occfilt_geo 'cellsize' method", {
  # Using cellsize method
  set.seed(1)
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("cellsize", factor = "3"),
    prj = crs(somevar)
  )
  expect_true(nrow(filtered) == 212)
})

test_that("occfilt_geo 'cellsize' method for different values", {
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("cellsize", factor = c(1, 4, 8)),
    prj = crs(somevar)
  )
  expect_equal(class(filtered), "list")
  expect_equal(length(filtered), 3)
})


test_that("occfilt_geo 'defined' method", {
  # Using defined method
  set.seed(1)
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("defined", d = "30"),
    prj = crs(somevar)
  )
  expect_true(nrow(filtered) == 78)
})

test_that("occfilt_geo 'defined' method for different values", {
  filtered <- occfilt_geo(
    data = spp1,
    x = "x",
    y = "y",
    env_layer = somevar,
    method = c("defined", factor = c(5, 10, 15)),
    prj = crs(somevar)
  )
  expect_equal(class(filtered), "list")
  expect_equal(length(filtered), 3)
})
