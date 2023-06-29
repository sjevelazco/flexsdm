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
  dplyr::slice_sample(n = 70)

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

test_that("test p_bpdp with continuous predictors and gam", {
  m_ <- fit_gam(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec"), k = -1
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and gau", {
  m_ <- fit_gau(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and glm", {
  m_ <- fit_glm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and gbm", {
  m_ <- fit_gbm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and max", {
  m_ <- fit_max(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and net", {
  m_ <- fit_net(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and raf", {
  m_ <- fit_raf(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})

test_that("test p_bpdp with continuous predictors and svm", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf), c("gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 0)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh", resolution = 5)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 25)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})


test_that("test p_bpdp with continuous and factor predictors and svm", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    predictors_f = "clusters",
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 8)

  # Partial dependence plot for training and projection condition found in a projection area
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50)
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})


test_that("training_boundaries = convexh and rectangle", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    predictors_f = "clusters",
    partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_bpdp(model = m_$model, training_data = abies2, training_boundaries = "convexh")
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(length(asdf$patches$plots), 8)

  # Partial dependence plot for training and projection condition and rectangle
  asdf <- p_bpdp(model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50, training_boundaries = "rectangle")
  expect_equal(class(asdf), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asdf$data), 2500)
})
