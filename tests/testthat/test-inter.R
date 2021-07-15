## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("original example 1- path difference 2-unable to find variable terra in example", {

    require(terra)
    require(dplyr)

    f <- system.file("external/suit_time_step.tif", package = "flexsdm")
    abma <- terra::rast(f)

    int <- inter(
      r1 = abma[[1]],
      r2 = abma[[2]],
      y1 = 2010,
      y2 = 2020,
      rastername = "Abies",
      dir_save = NULL
    )
    expect_true(class(int)=="SpatRaster")
    expect_equal(terra::nlyr(int), 11)
})

test_that("test save raster", {

  require(terra)
  require(dplyr)

  f <- system.file("external/suit_time_step.tif", package = "flexsdm")
  abma <- terra::rast(f)

  expect_message(int <- inter(
    r1 = abma[[1]],
    r2 = abma[[2]],
    y1 = 2010,
    y2 = 2020,
    rastername = "Abies",
    dir_save = tempdir()
  )
  )
})
