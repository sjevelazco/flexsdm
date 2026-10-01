test_that("sample_pseudoabs", {
  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  regions <- system.file("external/regions.tif", package = "flexsdm")
  regions <- terra::rast(regions)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp3") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # Pseudo-absences randomly sampled throughout study area
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = "random",
      rlayer = regions,
      maskval = NULL
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences k-means approach
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c(method = "kmeans", env = somevar),
      rlayer = regions,
      maskval = NULL
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences randomly sampled within a regions where a species occurs
  ## Regions where this species occurrs
  samp_here <- terra::extract(regions, single_spp[2:3])[, 2] %>%
    unique() %>%
    na.exclude()

  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = "random",
      rlayer = regions,
      maskval = samp_here
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences sampled with geographical constraint
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("geo_const", width = "30000"),
      rlayer = regions,
      maskval = samp_here
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences sampled with geo_const_kmeans constraint
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("geo_const_kmeans", width = "30000", env = somevar),
      rlayer = regions,
      maskval = samp_here
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences sampled with environmental constraint
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("env_const", env = somevar),
      rlayer = crop(regions, terra::ext(regions) - 33000),
      maskval = samp_here
    )

  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences sampled with environmental constraint and k-means
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("env_const_kmeans", env = somevar),
      rlayer = crop(regions, terra::ext(regions) - 33000),
      maskval = samp_here
    )

  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Pseudo-absences sampled with environmental and geographical constraint
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("geoenv_const", width = "50000", env = somevar),
      rlayer = crop(regions, terra::ext(regions) - 33000),
      maskval = samp_here
    )
  expect_equal(class(ps1)[1], "tbl_df")
  expect_equal(nrow(ps1), nrow(single_spp) * 10)
  rm(ps1)

  # Pseudo-absences sampled with environmental and geographical constraint and with k-mean
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("geoenv_const_kmeans", width = "50000", env = somevar),
      rlayer = crop(regions, terra::ext(regions) - 33000),
      maskval = samp_here
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  # Sampling pseudo-absence using a calibration area
  ca_ps1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("buffer", width = 50000),
    crs = crs(somevar)
  )

  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 50,
      method = "random",
      rlayer = regions,
      maskval = NULL,
      calibarea = ca_ps1
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)

  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 50,
      method = "random",
      rlayer = regions,
      maskval = samp_here,
      calibarea = ca_ps1
    )
  expect_equal(class(ps1)[1], "tbl_df")
  rm(ps1)
})


test_that("function misuse", {
  skip_on_cran()
  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  regions <- system.file("external/regions.tif", package = "flexsdm")
  regions <- terra::rast(regions)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp3") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # Pseudo-absences randomly sampled throughout study area
  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "env_const_XXXX",
    rlayer = regions,
    maskval = NULL
  ))

  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "kmeans",
    rlayer = regions,
    maskval = NULL
  ))

  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "env_const",
    rlayer = regions,
    maskval = NULL
  ))

  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "genv_const",
    rlayer = regions,
    maskval = NULL
  ))

  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "asdf_const_kmeans",
    rlayer = regions,
    maskval = NULL
  ))
})


test_that("sample_pseudoabs kmeans method works with maskval set", {
  skip_on_cran()
  # Regression test for issue #472: the kmeans branch masked rlayer using env
  # as the mask and assigned the result back into env (terra::mask(rlayer,
  # env) masks its FIRST argument), corrupting the environmental data handed
  # to kmeans with the categorical region raster instead. This only surfaces
  # when maskval is non-NULL (existing tests only covered maskval = NULL).
  set.seed(1)
  env <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 0, xmax = 20, ymin = 0, ymax = 20,
    nlyrs = 2
  )
  terra::values(env) <- cbind(rnorm(400), rnorm(400))
  names(env) <- c("v1", "v2")
  env[1] <- NA # keep at least one NA cell so kf()'s cell lookup is exercised normally

  rlayer <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 0, xmax = 20, ymin = 0, ymax = 20
  )
  terra::values(rlayer) <- rep(c(0, 1, 2), length.out = 400)
  levels(rlayer) <- data.frame(
    ID = c(0, 1, 2),
    category = c("water", "forest", "urban")
  )

  pts <- data.frame(x = runif(10, 0, 20), y = runif(10, 0, 20))

  ps1 <- sample_pseudoabs(
    data = pts,
    x = "x",
    y = "y",
    n = 10,
    method = c("kmeans", env = env),
    rlayer = rlayer,
    maskval = "forest"
  )

  expect_equal(class(ps1)[1], "tbl_df")
  expect_equal(nrow(ps1), 10)
  sampled_category <- suppressWarnings(
    terra::extract(rlayer, ps1[, c("x", "y")])$category
  )
  expect_true(all(sampled_category == "forest"))
})


test_that("sample_pseudoabs kmeans method works when env raster has no NA cells", {
  skip_on_cran()
  # Regression test for issue #472: kf() built cell ids from
  # names(km$cluster), but kmeans() receives a matrix produced via
  # dplyr::select(), which drops the data frame's row names, so
  # names(km$cluster) was always NULL. This previously only surfaced when
  # the environmental raster had no NA cells to filter (the common case is
  # masked because most real rasters have edge/ocean NAs).
  set.seed(1)
  env <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 0, xmax = 20, ymin = 0, ymax = 20,
    nlyrs = 2
  )
  terra::values(env) <- cbind(rnorm(400), rnorm(400)) # no NA cells

  rlayer <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 0, xmax = 20, ymin = 0, ymax = 20
  )
  terra::values(rlayer) <- 1

  pts <- data.frame(x = runif(10, 0, 20), y = runif(10, 0, 20))

  ps1 <- sample_pseudoabs(
    data = pts,
    x = "x",
    y = "y",
    n = 10,
    method = c("kmeans", env = env),
    rlayer = rlayer,
    maskval = NULL
  )

  expect_equal(class(ps1)[1], "tbl_df")
  expect_equal(nrow(ps1), 10)
})
