test_that("correct_colinvar Pearson", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  var <-
    correct_colinvar(
      env_layer = somevar,
      method = c("pearson", th = "0.8")
    )

  expect_equal(class(var), "list")
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(var$removed_variables, "CFP_3")
  expect_equal(nrow(var$correlation_table), 4)
  expect_true(all(names(var) %in% c("env_layer", "removed_variables", "correlation_table")))
})


test_that("correct_colinvar VIF", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  var <-
    correct_colinvar(env_layer = somevar, method = c("vif", th = "8"))

  expect_equal(class(var), "list")
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(var$removed_variables, "CFP_2")
  expect_equal(nrow(var$correlation_table), 3)
  expect_true(all(names(var) %in% c("env_layer", "removed_variables", "correlation_table")))
})


test_that("correct_colinvar PCA", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  var <-
    correct_colinvar(env_layer = somevar, method = "pca")

  expect_equal(length(var), 3)
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(nrow(var$coefficients), 4)
  expect_equal(nrow(var$cumulative_variance), 4)
  expect_true(all(names(var) %in% c("env_layer", "coefficients", "cumulative_variance")))
})

test_that("correct_colinvar FA", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  var <- correct_colinvar(env_layer = somevar, method = c("fa"))

  expect_equal(length(var), 4)
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(length(var$removed_variables), 3)
  expect_equal(nrow(var$correlation_table), 4)
  expect_equal(ncol(var$variable_loading), 2)
  expect_true(all(names(var) %in% c("env_layer", "removed_variables", "correlation_table", "variable_loading")))
})
