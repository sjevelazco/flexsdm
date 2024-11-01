library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field = "clusters")
somevar <- c(somevar, cat)

abies2 <- abies %>%
  select(x, y, pr_ab)

abies2 <- sdm_extract(abies2,
  x = "x",
  y = "y",
  env_layer = somevar
)
abies2 <- part_random(abies2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)

svm_t1 <- fit_svm(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  predictors_f = "clusters",
  partition = ".part",
  thr = c("max_sens_spec")
)

test_that("test data_bpdp with factor", {
  df <- data_bpdp(
    model = svm_t1$model,
    predictors = c("cwd", "clusters"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(names(df), c("pspdata", "training_boundaries"))
  expect_equal(nrow(df$pspdata), 1000)
})

test_that("test data_bpdp", {
  df <- data_bpdp(
    model = svm_t1$model,
    predictors = c("aet", "cwd"),
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(names(df), c("pspdata", "training_boundaries"))
  expect_equal(length(df), 2)

  df <- data_bpdp(
    model = svm_t1$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(nrow(df$pspdata), 2500)

  df <- data_bpdp(
    model = svm_t1$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = NULL,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(nrow(df$pspdata), 2500)
  expect_true(is.null(df$training_boundaries))

  df <- data_bpdp(
    model = svm_t1$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    training_boundaries = "convexh",
    projection_data = NULL,
    training_data = abies2,
    clamping = FALSE
  )

  expect_true(nrow(df$training_boundaries) == 13)
})


library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field = "clusters")
somevar <- c(somevar, cat)

data(abies)

set.seed(123)
abies2 <- abies %>%
  select(x, y, pr_ab) %>%
  group_by(pr_ab) %>%
  dplyr::slice_sample(n = 50)

set.seed(210)
abies2 <- sdm_extract(abies2,
  x = "x",
  y = "y",
  env_layer = somevar
) %>%
  part_random(
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )



test_that("test psp with gam", {
  m_ <- fit_gam(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    training_boundaries = NULL,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})

test_that("test psp with gau", {
  m_ <- fit_gau(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "rectangle",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 4)
})

test_that("test psp with gbm", {
  m_ <- fit_gbm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "rectangle",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 4)
})

test_that("test psp with glm", {
  m_ <- fit_glm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})

test_that("test psp with max", {
  m_ <- fit_max(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})

test_that("test psp with net", {
  m_ <- fit_net(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "rectangle",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 4)
})

test_that("test psp with raf", {
  m_ <- fit_raf(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})

test_that("test psp with svm", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})


test_that("test psp with factors", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "cwd"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pspdata), 2500)
  expect_equal(length(df), 2)
  expect_true(is.null(df$training_boundaries))

  # with convex hull
  df <- data_bpdp(
    model = m_$model,
    predictors = c("aet", "tmx"),
    resolution = 50,
    projection_data = somevar,
    training_data = abies2,
    training_boundaries = "convexh",
    clamping = FALSE
  )
  expect_true(nrow(df$training_boundaries) == 8)
})
