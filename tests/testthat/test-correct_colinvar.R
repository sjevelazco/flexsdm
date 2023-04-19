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
      method = c("pearson", th = "0.5")
    )

  expect_equal(class(var), "list")
  expect_equal(length(var$cor_variables), 4)
  expect_equal(nrow(var$cor_table), 4)
  expect_true(all(names(var) %in% c("cor_table", "cor_variables")))

  # test default correlation threshold
  var <- correct_colinvar(env_layer = somevar,
                   method = c("pearson"))
  expect_equal(length(var$cor_variables), 4)

  # test with maxcell
  var <- correct_colinvar(env_layer = somevar,
                          method = c("pearson"), maxcell = 5000)
  expect_equal(length(var$cor_variables), 4)

  # To high th, no variable is expected to be removed
  var <- correct_colinvar(env_layer = somevar,
                          method = c("pearson", th = "0.9"))
  expect_equal(var$cor_variables, "No pair of variables reached the specified correlation threshold.")
})


test_that("correct_colinvar VIF", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  CPF_5 <- somevar[[2]]
  names(CPF_5) <- 'CPF_5'
  somevar <- terra::rast(list(somevar, somevar[[2]]))
  names(somevar)[5] <- "CPF_5"
  var <-
    correct_colinvar(env_layer = somevar, method = c("vif", th = "8"))

  expect_equal(class(var), "list")
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(var$removed_variables, "CFP_2")
  expect_true(all(names(var) %in% c("env_layer", "removed_variables", "vif_table")))

  # test with maxcell
  var <-
    correct_colinvar(env_layer = somevar, method = c("vif", th = "8"), maxcell = 50000)

  expect_equal(class(var), "list")
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(var$removed_variables, "CFP_2")
  expect_true(all(names(var) %in% c("env_layer", "removed_variables", "vif_table")))
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

   # Test with maxcell
  var <-
    correct_colinvar(env_layer = somevar, method = "pca", maxcell = 5000)

  expect_equal(length(var), 3)
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(nrow(var$coefficients), 4)
  expect_equal(nrow(var$cumulative_variance), 4)
  expect_true(all(names(var) %in% c("env_layer", "coefficients", "cumulative_variance")))
})

test_that("correct_colinvar PCA with projections", {
  require(terra)
  require(dplyr)
  dir_sc <- file.path(tempdir(), "projections")
  dir.create(dir_sc)
  dir_sc <- file.path(dir_sc, c('scenario_1', 'scenario_2'))
  sapply(dir_sc, dir.create)

  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  terra::writeRaster(somevar, file.path(dir_sc[1], "somevar.tif"), overwrite=TRUE)
  terra::writeRaster(somevar, file.path(dir_sc[2], "somevar.tif"), overwrite=TRUE)

  # Perform pearson collinearity control
  var <-
    correct_colinvar(env_layer = somevar, method = "pca", proj = dirname(dir_sc[1]))

  expect_true(is.character(var$proj))
  unlink(dirname(dir_sc[1]), recursive = TRUE)
  unlink(gsub("projections", "Projection_PCA", dirname(dir_sc[1])), recursive = TRUE)
})

test_that("correct_colinvar FA", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Perform pearson collinearity control
  var <- correct_colinvar(env_layer = somevar, method = c("fa"))

  expect_equal(length(var), 5)
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(length(var$removed_variables), 3)
  expect_equal(nrow(var$loadings), 4)
  expect_true(all(names(var) %in% c("env_layer", "number_factors", "removed_variables", "uniqueness", "loadings")))

  # Test with maxcell
  var <- correct_colinvar(env_layer = somevar, method = c("fa"), maxcell = 5000)

  expect_equal(length(var), 5)
  expect_equal(class(var$env_layer)[1], "SpatRaster")
  expect_equal(length(var$removed_variables), 3)
  expect_equal(nrow(var$loadings), 4)
  expect_true(all(names(var) %in% c("env_layer", "number_factors", "removed_variables", "uniqueness", "loadings")))
})

test_that("misuse of argument", {
  require(terra)
  require(dplyr)
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")

  # Perform pearson collinearity control
  expect_error(var <- correct_colinvar(env_layer = somevar, method = c("faasdf")))
})
