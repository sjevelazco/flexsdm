test_that("get_block", {
  require(dplyr)
  require(terra)
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Lest practice with a single species
  single_spp <- spp %>% dplyr::filter(species == "sp3")

  part <- part_sblock(
    env_layer = somevar,
    data = single_spp,
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
  expect_equal(class(grid_env)[[1]], "SpatRaster")
})
