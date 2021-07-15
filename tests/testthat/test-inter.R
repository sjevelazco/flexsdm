## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("original example 1- path difference 2-unable to find variable terra in example", {

    require(terra)
    require(dplyr)

    f <- system.file("/inst/external/suit_time_step.tif", package = "flexsdm")
    abma <- terra::rast(f)
    plot(abma)

    expect_error(int <- inter(
      r1 = abma[[1]],
      r2 = abma[[2]],
      y1 = 2010,
      y2 = 2020,
      rastername = "Abies",
      dir_save = NULL,
      n_cores = 1
    ))
    int
})

