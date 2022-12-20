library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field="clusters")
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

test_that("test data_pdp with factor", {

  df <- data_pdp(
    model = svm_t1$model,
    predictors = c("clusters"),
    resolution = 100,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(names(df), c("pdpdata", "resid"))
  expect_equal(nrow(df$pdpdata), 14)
})

test_that("test data_pdp", {

  df <- data_pdp(
    model = svm_t1$model,
    predictors = c("aet"),
    resolution = 100,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(names(df), c("pdpdata", "resid"))
  expect_equal(length(df), 2)

  df <- data_pdp(
    model = svm_t1$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(nrow(df$pdpdata), 50)

  df <- data_pdp(
    model = svm_t1$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = NULL,
    training_data = abies2,
    clamping = FALSE
  )

  expect_equal(nrow(df$pdpdata), 50)

  df <- data_pdp(
    model = svm_t1$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = NULL,
    training_data = abies2,
    clamping = FALSE
  )

  expect_true(is.null(df$resid))
})


library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field="clusters")
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



test_that("test pdp with gam", {
  m_ <- fit_gam(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with gau", {
  m_ <- fit_gau(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with gbm", {
  m_ <- fit_gbm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with glm", {
  m_ <- fit_glm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with max", {
  m_ <- fit_max(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with net", {
  m_ <- fit_net(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)
  expect_error(data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = NULL,
    clamping = FALSE
  ))

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with raf", {
  m_ <- fit_raf(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})

test_that("test pdp with svm", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})


test_that("test pdp with factors", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), partition = ".part", thr = c("max_sens_spec")
  )
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = TRUE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_equal(nrow(df$pdpdata), 50)
  expect_equal(length(df), 2)

  # without residuals
  df <- data_pdp(
    model = m_$model,
    predictors = c("aet"),
    resolution = 50,
    resid = FALSE,
    projection_data = somevar,
    training_data = abies2,
    clamping = FALSE
  )
  expect_true(is.null(df$resid))
})




