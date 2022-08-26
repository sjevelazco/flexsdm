test_that("extra_truncate", {
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
  ca <-
    calib_area(
      sp,
      x = "x",
      y = "y",
      method = c("buffer", width = 30000),
      crs = crs(somevar)
    )

  # Get environmental condition of calibration area
  somevar_ca <- somevar %>%
    crop(., ca) %>%
    mask(., ca)


  xp <-
    extra_eval(
      training_data = somevar_ca,
      projection_data = somevar,
      n_cores = 2,
      aggreg_factor = 5
    )
  xp2 <- extra_truncate(
    suit = somevar$CFP_1,
    extra = xp,
    threshold = 50
  )
  expect_equal(class(xp2[[1]])[1], "SpatRaster")
})
