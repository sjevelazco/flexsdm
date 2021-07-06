test_that("multiplication works", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <-
    system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  set.seed(1)
  some_sp <- spp %>%
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
  mglm <- fit_glm(
    data = some_sp,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
    poly = 2
  )
  mraf <- fit_raf(
    data = some_sp,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part",
  )
  mgbm <- fit_gbm(
    data = some_sp,
    response = "pr_ab",
    predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
    partition = ".part"
  )

  # Fit and ensemble
  mensemble <- fit_ensemble(
    models = list(mglm, mraf, mgbm),
    ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
    thr = NULL,
    thr_model = "max_sens_spec",
    metric = "TSS"
  )

  expect_equal(class(mensemble), 'list')
  expect_equal(length(mensemble), 4)
  expect_equal(
    unique(mensemble$performance$model),
    c("mean", "meanw", "meansup", "meanthr", "median")
  )


  # Expect an error
  expect_error(fit_ensemble(
    models = list(mglm, mraf, mgbm),
    ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
    thr = NULL,
    thr_model = NULL,
    metric = NULL
  ))
})
