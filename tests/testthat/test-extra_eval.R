test_that("extra_eval one and two cores", {
  require(dplyr)
  require(terra)

  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  spp$species %>% unique()
  sp <- spp %>%
    dplyr::filter(species == "sp3", pr_ab == 1) %>%
    dplyr::select(x, y)

  # Accessible area
  ca <- calib_area(sp, x = "x", y = "y", method = c("buffer", width = 30000), crs=crs(somevar))

  # Get environmental condition of calibration area
  somevar_ca <- somevar %>%
    crop(., ca) %>%
    mask(., ca)


  xp <-
    extra_eval(
      env_calib = somevar_ca,
      env_proj = somevar,
      n_cores = 2,
      aggreg_factor = 3
    )
  expect_equal(class(xp)[1], "SpatRaster")

  somevar_ca <- terra::aggregate(somevar_ca, 3)
  somevar <- terra::aggregate(somevar, 3)

  xp2 <-
    extra_eval(
      env_calib = somevar_ca,
      env_proj = somevar,
      n_cores = 1,
      aggreg_factor = 1
    )
  expect_equal(class(xp2)[1], "SpatRaster")
})


test_that("extra_eval wrong use", {
  require(dplyr)
  require(terra)

  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  spp$species %>% unique()
  sp <- spp %>%
    dplyr::filter(species == "sp3", pr_ab == 1) %>%
    dplyr::select(x, y)

  # Accessible area
  ca <- calib_area(sp, x = "x", y = "y", method = c("buffer", width = 30000), crs=crs(somevar))

  # Get environmental condition of calibration area
  somevar_ca <- somevar %>%
    crop(., ca) %>%
    mask(., ca)

  names(somevar_ca) <- tolower(names(somevar_ca))

  expect_error(extra_eval(
    env_calib = somevar_ca,
    env_proj = somevar,
    n_cores = 1,
    aggreg_factor = 3
  ))
})
