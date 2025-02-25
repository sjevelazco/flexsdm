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
      method = c(method = 'kmeans', env = somevar),
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

  # Pseudo-absences sampled with environmental and geographical constraint
  ps1 <-
    sample_pseudoabs(
      data = single_spp,
      x = "x",
      y = "y",
      n = nrow(single_spp) * 10,
      method = c("geo_env_const", width = "50000", env = somevar),
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
      method = c("geo_env_km_const", width = "50000", env = somevar),
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
    method = "geo_env_const",
    rlayer = regions,
    maskval = NULL
  ))

  expect_error(sample_pseudoabs(
    data = single_spp,
    x = "x",
    y = "y",
    n = nrow(single_spp) * 10,
    method = "geo_env_km_const",
    rlayer = regions,
    maskval = NULL
  ))
})
