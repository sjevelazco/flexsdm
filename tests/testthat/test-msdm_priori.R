test_that("msdm_priori", {
  require(dplyr)
  require(terra)

  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # It will be select the presences of a species
  occ <- spp %>%
    dplyr::filter(species == "sp3", pr_ab == 1)

  # It will select a raster layer to be used as a basic raster
  a_variable <- somevar[[1]]

  ### xy method
  m_xy <- msdm_priori(
    data = occ,
    x = "x",
    y = "y",
    method = "xy",
    env_layer = a_variable
  )

  expect_s4_class(m_xy, "SpatRaster")

  ### min method
  m_min <- msdm_priori(
    data = occ,
    x = "x",
    y = "y",
    method = "min",
    env_layer = a_variable
  )
  expect_s4_class(m_min, "SpatRaster")

  ### cml method
  m_cml <- msdm_priori(
    data = occ,
    x = "x",
    y = "y",
    method = "cml",
    env_layer = a_variable
  )
  expect_s4_class(m_cml, "SpatRaster")

  ### ker method
  m_ker <- msdm_priori(
    data = occ,
    x = "x",
    y = "y",
    method = "ker",
    env_layer = a_variable
  )
  expect_s4_class(m_ker, "SpatRaster")

  # ### Object different to SpatRaster
  # require(raster)
  # m_ker <- msdm_priori(
  #   data = occ,
  #   x = "x",
  #   y = "y",
  #   method = "ker",
  #   env_layer = a_variable
  # )
  # expect_s4_class(m_ker, "SpatRaster")
})


test_that("function misuse", {
  require(dplyr)
  require(terra)

  data("spp")
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # It will be select the presences of a species
  occ <- spp %>%
    dplyr::filter(species == "sp3", pr_ab == 1)

  # It will select a raster layer to be used as a basic raster
  a_variable <- somevar[[1]]

  expect_error(msdm_priori(
    data = occ,
    x = "x",
    y = "y",
    method = "xy",
    env_layer = NULL
  ))

  expect_error(msdm_priori(
    data = occ,
    # x = "x",
    y = "y",
    method = "xy",
    env_layer = a_variable
  ))

  expect_error(msdm_priori(
    data = occ,
    y = "y",
    method = "xy",
    env_layer = a_variable
  ))
})
