test_that("conventional partition situtation", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Lest practice with a single species
  single_spp <- spp %>% dplyr::filter(species == "sp3")
  # plot(somevar$CFP_1)
  # single_spp[c('x', 'y')] %>% points

  single_spp$x[which.max(single_spp$x)] <- -150540.8
  single_spp$x[1:5] <- NA
  part <- part_sblock(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 10,
    max_res_mult = 500,
    num_grids = 30,
    n_part = 2,
    min_occ = 5,
    prop = 0.5
  )
  expect_equal(class(part), "list")
  expect_equal(length(part), 3)

  single_spp
  expect_message(part_sblock(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 10,
    max_res_mult = 500,
    num_grids = 30,
    n_part = 2,
    prop = 0.5
  ))
})


test_that("only with presences", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Lest practice with a single species
  single_spp <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
  # plot(somevar$CFP_1)
  # single_spp[c('x', 'y')] %>% points

  single_spp$x[which.max(single_spp$x)] <- -150540.8
  part <- part_sblock(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 10,
    max_res_mult = 500,
    num_grids = 30,
    n_part = 2,
    prop = 0.5
  )
  expect_equal(class(part), "list")
  expect_equal(length(part), 3)
})


test_that("tese some errors", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Lest practice with a single species
  single_spp <- spp %>% dplyr::filter(species == "sp3", pr_ab == 1)
  single_spp$pr_ab[1:3] <- 1:3

  expect_error(part <- part_sblock(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 10,
    max_res_mult = 500,
    num_grids = 30,
    n_part = 2,
    min_occ = 10,
    prop = 0.5
  ))
})
