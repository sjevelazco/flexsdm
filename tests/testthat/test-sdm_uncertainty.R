# Resampling Methods
require(dplyr)
require(terra)

data("spp")
somevar <- system.file("external/somevar.tif", package = "flexsdm") %>%
  terra::rast()

# Extract data
some_sp <- spp %>%
  filter(species == "sp3")

some_sp <-
  sdm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = somevar
  )

data("backg")
backg <-
  sdm_extract(
    data = backg[1:3],
    x = "x",
    y = "y",
    env_layer = somevar
  )
names(backg)[1] <- "pres_abs"

some_sp <- rename(some_sp, "pres_abs" = "pr_ab")

# Partition
some_sp <- part_random(
  data = some_sp,
  pr_ab = "pres_abs",
  method = c(method = "kfold", folds = 3)
)
backg <- part_random(
  data = backg,
  pr_ab = "pres_abs",
  method = c(method = "kfold", folds = 3)
)

# Resample to tests speed up
somevar <- terra::aggregate(somevar, 6)

test_that("test with GAM", {
  m <- fit_gam(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part", k = 2
  )
  unc <- sdm_uncertainty(
    model = m,
    training_data = some_sp,
    response = "pres_abs",
    projection_data = somevar,
    iteration = 5,
    n_cores = 3
  )

  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with GLM", {
  m <- fit_glm(
  data = some_sp,
  response = "pres_abs",
  predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
  partition = ".part",
  poly = 2,
  inter_order = 2
)
  unc <- sdm_uncertainty(
    model = m,
    training_data = some_sp,
    response = "pres_abs",
    projection_data = somevar,
    iteration = 5,
    n_cores = 3
  )

  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with GBM", {
  m <- tune_gbm(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    grid = expand.grid(
    n.trees = c(50, 100),
    shrinkage = c(1),
    n.minobsinnode = c(5, 9)
  ),
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 1
  )
  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with GAU", {
  m <- fit_gau(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    background = NULL,
    partition = ".part"
  )
  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    background = NULL,
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with NET", {
  m <- tune_net(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    grid = expand.grid(
    size = c(2:3),
    decay = c(1, 3)
  ),
  n_cores = 2,
  thr = "equal_sens_spec"
  )
  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    background = NULL,
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with RAF", {
  m <- tune_raf(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    grid = expand.grid(
    mtry = seq(1, 7, 1),
    ntree = c(400, 600, 800)
  ),
  n_cores = 2,
  thr = "equal_sens_spec"
  )
  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    background = NULL,
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with SVM", {
  m <- tune_svm(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    grid = expand.grid(
    C = c(2, 8, 20),
    sigma = c(0.01, 0.1, 0.4)
  ),
  n_cores = 2,
  thr = "equal_sens_spec"
  )
  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    background = NULL,
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with MAX", {
  m <- tune_max(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    background = backg,
    grid = expand.grid(
    regmult = seq(0.1, 3, 0.5),
    classes = c("lqpht")
  ),
  n_cores = 2,
  thr = "equal_sens_spec"
  )

  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    background = backg,
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

test_that("test with DOM", {
  m <- fit_dom(
    data = some_sp,
    response = "pres_abs",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part"
  )

  unc <- sdm_uncertainty(
    models = m,
    training_data = some_sp,
    response = "pres_abs",
    projection_data = somevar,
    iteration = 5,
    n_cores = 2
  )
  expect_equal(class(unc)[[1]], "SpatRaster")
  expect_true(terra::global(unc, max, na.rm = TRUE)[1, 1] > 0)
})

