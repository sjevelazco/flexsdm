test_that("msdm_posteriori", {
  require(dplyr)
  require(terra)

  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)


  # It will prepared data for modeling a species
  set.seed(10)
  occ <- spp %>%
    dplyr::filter(species == "sp2") %>% # filter a species
    sdm_extract(
      data = ., x = "x", y = "y",
      env_layer = somevar, filter_na = TRUE
    ) %>% # extrac variables values
    part_random(.,
                pr_ab = "pr_ab",
                method = c(method = "kfold", folds = 10)
    ) # add columns with partition

  # Lets fit and predict a model
  m_glm <- fit_glm(
    data = occ,
    response = "pr_ab",
    predictors = names(somevar),
    partition = ".part",
    thr = "equal_sens_spec",
  )

  # Lets predict this model
  m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)

  ### bmcp method
  m_bmcp <- msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = 30000
  )
  expect_s4_class(m_bmcp, "SpatRaster")

  ### mcp method
  m_mcp <- msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "mcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = NULL
  )

  expect_s4_class(m_mcp, "SpatRaster")


  ### pres method
  m_pres <- msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "pres",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = NULL
  )

  expect_s4_class(m_pres, "SpatRaster")

  ### lq method
  m_lq <- msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "lq",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = NULL
  )
  expect_s4_class(m_lq, "SpatRaster")


  ### obr method
  m_obr <- msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "obr",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = NULL
  )
  expect_s4_class(m_obr, "SpatRaster")
})

test_that("use of other object class", {
  require(dplyr)
  require(terra)

  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)


  # It will prepared data for modeling a species
  set.seed(10)
  occ <- spp %>%
    dplyr::filter(species == "sp2") %>% # filter a species
    sdm_extract(
      data = ., x = "x", y = "y",
      env_layer = somevar, filter_na = TRUE
    ) %>% # extrac variables values
    part_random(.,
                pr_ab = "pr_ab",
                method = c(method = "kfold", folds = 10)
    ) # add columns with partition

  # Lets fit and predict a model
  m_glm <- fit_glm(
    data = occ,
    response = "pr_ab",
    predictors = names(somevar),
    partition = ".part",
    thr = "equal_sens_spec",
  )

  # Lets predict this model
  m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)

  # No use of buffer with bmcp
  expect_s4_class(msdm_posteriori(
    records = data.frame(occ),
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = raster::raster(m_pred[[1]]),
    thr = "equal_sens_spec",
    buffer = 30000
  ), "SpatRaster")
})

test_that("missuse of function", {
  require(dplyr)
  require(terra)

  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)


  # It will prepared data for modeling a species
  set.seed(10)
  occ <- spp %>%
    dplyr::filter(species == "sp2") %>% # filter a species
    sdm_extract(
      data = ., x = "x", y = "y",
      env_layer = somevar, filter_na = TRUE
    ) %>% # extrac variables values
    part_random(.,
                pr_ab = "pr_ab",
                method = c(method = "kfold", folds = 10)
    ) # add columns with partition

  # Lets fit and predict a model
  m_glm <- fit_glm(
    data = occ,
    response = "pr_ab",
    predictors = names(somevar),
    partition = ".part",
    thr = "equal_sens_spec",
  )

  # Lets predict this model
  m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)

  # No use of buffer with bmcp
  expect_error(msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec",
    buffer = NULL
  ))

  # miss use of threshold
  expect_error(msdm_posteriori(
    records = occ,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec_XXXX",
    buffer = 30000
  ))

  # miss use of x y
  expect_error(msdm_posteriori(
    records = occ,
    x = "xXXXX",
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec_XXXX",
    buffer = 30000
  ))

  # miss use of x y
  expect_error(msdm_posteriori(
    records = occ,
    x = NA,
    y = "y",
    pr_ab = "pr_ab",
    method = "bmcp",
    cont_suit = m_pred[[1]],
    thr = "equal_sens_spec_XXXX",
    buffer = 30000
  ))
})
