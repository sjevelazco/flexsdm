test_that("test for fit_ function family", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  set.seed(1)
  spp_ <- spp %>%
    dplyr::filter(species == "sp2") %>%
    sdm_extract(
      data = .,
      x = 'x',
      y = 'y',
      env_layer = somevar,
      variables = names(somevar),
      filter_na = TRUE
    ) %>%
    part_random(
      data = .,
      pr_ab = 'pr_ab',
      method = c(method = "kfold", folds = 3)
    )


  # gam
  gam <- fit_gam(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = gam,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )

  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # gau
  gau <- fit_gau(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = gau,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # gbm
  gbm <- fit_gbm(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = gbm,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # glm
  glm <- fit_glm(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = glm,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # max
  max <- fit_max(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = max,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE,
    clamp = TRUE,
    pred_type = "cloglog"
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # net
  net <- fit_net(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec"),
    size = 1
  )

  p <- sdm_predict(
    models = net,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # raf
  raf <- fit_raf(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = raf,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # svm
  svm <- fit_svm(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  p <- sdm_predict(
    models = svm,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(class(p[[1]]) == "SpatRaster")
  rm(p)

  # Predict list of individual models
  p <- sdm_predict(
    models = list(svm, raf),
    pred = somevar,
    thr = NULL,
    con_thr = FALSE
  )
  expect_true(length(p) == 2)
  expect_equal(names(p) , c('svm', 'raf'))
  rm(p)
})


test_that("test for ensemble, mask, and suit. values above threshold", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  set.seed(1)
  spp_ <- spp %>%
    dplyr::filter(species == "sp2") %>%
    sdm_extract(
      data = .,
      x = 'x',
      y = 'y',
      env_layer = somevar,
      variables = names(somevar),
      filter_na = TRUE
    ) %>%
    part_random(
      data = .,
      pr_ab = 'pr_ab',
      method = c(method = "kfold", folds = 3)
    )

  ca <- calib_area(data = spp_, 'x', 'y', method = 'mcp')

  # gau
  gau <- fit_gau(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part"
  )

  # gbm
  gbm <- fit_gbm(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part"
  )

  # glm
  glm <- fit_glm(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    predictors_f = NULL,
    partition = ".part"
  )

  enm <-
    fit_ensemble(models = list(gau, gbm, glm),
                 ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
                 metric="TSS",
                 thr_model="equal_sens_spec"
                 )

  # Test predict ensemble and with predict_area
  p <- sdm_predict(
    models = enm,
    pred = somevar,
    thr = NULL,
    con_thr = FALSE,
    predict_area = ca
  )

  expect_true(class(p[[1]]) == "SpatRaster")
  expect_true(length(p) == 5)
  expect_false(terra::ext(p[[1]]) == ext(somevar))

  # Test predict ensemble and with predict_area and con_thr = TRUE
  p <- sdm_predict(
    models = enm,
    pred = somevar,
    thr = "max_sens_spec",
    predict_area = ca,
    con_thr = TRUE
  )

  expect_true(nrow(unique(p[[1]][[2]])) > 100)
  expect_false(terra::ext(p[[1]]) == ext(somevar))
})



test_that("test for all threshold", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  set.seed(1)
  spp_ <- spp %>%
    dplyr::filter(species == "sp2") %>%
    sdm_extract(
      data = .,
      x = 'x',
      y = 'y',
      env_layer = somevar,
      variables = names(somevar),
      filter_na = TRUE
    ) %>%
    part_random(
      data = .,
      pr_ab = 'pr_ab',
      method = c(method = "kfold", folds = 3)
    )

  ca <- calib_area(data = spp_, 'x', 'y', method = 'mcp')

  gam <- fit_gam(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part"
  )

  p <- sdm_predict(
    models = gam,
    pred = somevar,
    thr = "all",
    predict_area = ca,
    con_thr = TRUE
  )

  expect_equal(terra::nlyr(p[[1]]), 8)
})


test_that("test for prdicting ensemble of small models", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  set.seed(1)
  spp_ <- spp %>%
    dplyr::filter(species == "sp2") %>%
    sdm_extract(
      data = .,
      x = 'x',
      y = 'y',
      env_layer = somevar,
      variables = names(somevar),
      filter_na = TRUE
    ) %>%
    part_random(
      data = .,
      pr_ab = 'pr_ab',
      method = c(method = "kfold", folds = 3)
    )

  ca <- calib_area(data = spp_, 'x', 'y', method = 'mcp')

  gam <- esm_gam(
    data = spp_,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part"
  )

  p <- sdm_predict(
    models = gam,
    pred = somevar,
    thr = NULL,
    predict_area = ca,
    con_thr = TRUE
  )

  expect_equal(terra::nlyr(p[[1]]), 1)
})
