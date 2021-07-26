test_that("test occfilt_env", {
  require(terra)
  require(dplyr)

  # Environmental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp1 <- spp %>% dplyr::filter(species == "sp2", pr_ab == 1)
  spp1$idd <- 1:nrow(spp1)
  spp1$y[1] <- -300000

  # 5 bins conventional without use of cores
  filtered_1 <- occfilt_env(
    data = spp1,
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar,
    nbins = 5,
    cores = 1
  )

  expect_true(nrow(filtered_1) < nrow(spp1))

  #
  #   # 5 bins conventional without use of cores
  #   filtered_2 <- occfilt_env(
  #     data = spp1,
  #     x = "x",
  #     y = "y",
  #     id = "idd",
  #     env_layer = somevar,
  #     nbins = 5,
  #     cores = 3
  #   )
  #
  #   expect_true(nrow(filtered_1)==nrow(filtered_2))
})
