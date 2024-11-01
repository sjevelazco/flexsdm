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

test_that("test occfilt_env", {
  # 5 bins conventional without use of cores
  filtered_1 <- occfilt_env(
    data = spp1,
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar,
    nbins = 5
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
  #     nbins = 5
  #    )
  #
  #   expect_true(nrow(filtered_1)==nrow(filtered_2))
})

test_that("test occfilt_env with different values", {
  filtered_1 <- occfilt_env(
    data = spp1,
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar,
    nbins = c(5, 7, 9)
  )

  expect_equal(class(filtered_1), "list")
  expect_equal(length(filtered_1), 3)
  expect_true(all(sapply(filtered_1, nrow) %in% c("16", "20", "21")))
})

test_that("test occfilt_env with factor as variables", {
  somevar <- system.file("external/somevar.tif", package = "flexsdm") %>% terra::rast()
  somevarf <- system.file("external/regions.tif", package = "flexsdm") %>% terra::rast()
  somevar <- c(somevar, somevarf)
  filtered_1 <- occfilt_env(
    data = spp1,
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar,
    nbins = 5
  )

  expect_true(nrow(filtered_1) == 16)
})
