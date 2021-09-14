test_that("sample_background random method", {
  require(dplyr)
  data(spp)
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Lest practice with a single species
  spp_pa <- spp %>% dplyr::filter(species == "sp3")

  #
  part <- part_sblock(
    env_layer = somevar,
    data = spp_pa,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 100,
    max_res_mult = 500,
    num_grids = 10,
    min_occ = 2,
    n_part = 2
  )
  grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
  spp_p <- spp_pa %>% dplyr::filter(pr_ab == 1)


  # Sample background points throughout study area with random method
  bg <-
    sample_background(
      data = spp_p,
      x = "x",
      y = "y",
      n = 1000,
      method = "random",
      rlayer = grid_env
    )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  # Sample random background points constrained to a region with a give set of values
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 1000,
    method = "random",
    rlayer = grid_env,
    maskval = 1
  )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 1000,
    method = "random",
    rlayer = grid_env,
    maskval = 2
  )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)


  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 1000,
    method = "random",
    rlayer = grid_env,
    maskval = c(1, 2)
  )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)


  # Sample random background within a calibration area and constrained to a region
  ca_ps1 <- calib_area(
    data = spp_pa,
    x = "x",
    y = "y",
    method = c("buffer", width = 50000),
  )
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 1000,
    method = "random",
    rlayer = grid_env,
    maskval = 1,
    calibarea = ca_ps1
  )
  bg <- expect_equal(class(bg)[1], "tbl_df")
  rm(bg)
})

test_that("sample_background Thickening method", {
  require(dplyr)
  data(spp)
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Lest practice with a single species
  spp_pa <- spp %>% dplyr::filter(species == "sp3")

  #
  part <- part_sblock(
    env_layer = somevar,
    data = spp_pa,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_res_mult = 100,
    max_res_mult = 500,
    num_grids = 10,
    min_occ = 2,
    n_part = 2
  )
  grid_env <- get_block(env_layer = somevar, best_grid = part$grid)

  spp_p <- spp_pa %>% dplyr::filter(pr_ab == 1)

  # Thickening background without containing them
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 5000,
    method = "thickening",
    rlayer = grid_env,
  )

  bg <- expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  # Thickening background without using a given buffer width
  spp_p # presences database of a species
  grid_env # The reaster layer used for sampling background
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 5000,
    method = c("thickening", width = 150000),
    rlayer = grid_env
  )

  bg <- expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  # Sample thickening background within a calibration area and constrained to a region
  ca_ps1 <- calib_area(
    data = spp_pa,
    x = "x",
    y = "y",
    method = c("buffer", width = 50000),
  )

  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 3000,
    method = "thickening",
    rlayer = grid_env,
    maskval = 2,
    calibarea = ca_ps1
  )
  bg <- expect_equal(class(bg)[1], "tbl_df")
  rm(bg)
})


test_that("sample_background biased method", {
  require(dplyr)
  require(terra)
  data(spp)

  # Lets select the presences of a species
  spp_p <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

  # Raster layer with density of poinst to obtain a biased sampling background
  occ_density <- system.file("external/occ_density.tif", package = "flexsdm")
  occ_density <- terra::rast(occ_density)

  # A layer with region used to contraint background
  regions <- system.file("external/regions.tif", package = "flexsdm")
  regions <- terra::rast(regions)


  # Biased background points
  spp_p # presences database of a species
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 3000,
    method = "biased",
    rlayer = regions,
    rbias = occ_density
  )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  # Biased background points constrained in a region
  bg <- sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 500,
    method = "biased",
    rlayer = regions,
    rbias = occ_density,
    maskval = c(1, 2)
  )
  expect_equal(class(bg)[1], "tbl_df")
  rm(bg)

  # Exeed number of points
  expect_message(sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 500000,
    method = "biased",
    rlayer = raster::raster(regions),
    rbias = raster::raster(occ_density),
    maskval = c(1, 2)
  ))
})

test_that("sample_background misuse of argument", {
  require(dplyr)
  require(terra)
  data(spp)

  # Lets select the presences of a species
  spp_p <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

  # Raster layer with density of poinst to obtain a biased sampling background
  occ_density <- system.file("external/occ_density.tif", package = "flexsdm")
  occ_density <- terra::rast(occ_density)

  # A layer with region used to contraint background
  regions <- system.file("external/regions.tif", package = "flexsdm")
  regions <- terra::rast(regions)


  # Biased background points
  expect_error(sample_background(
    data = spp_p,
    x = "x",
    y = "y",
    n = 3000,
    method = "biasdfsed",
    rlayer = regions,
    rbias = occ_density
  ))

  expect_error(
    sample_background(
      data = spp_p,
      x = "x",
      y = "y",
      n = 3000,
      method = "biased",
      rlayer = regions
      # rbias = occ_density
    ))
})
