test_that("part_sband lon", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Let's practice with two longitudinal partition with presences and absences
  single_spp <- spp %>% dplyr::filter(species == "sp1")
  part_1 <- part_sband(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    type = "lon",
    min_bands = 2,
    max_bands = 20,
    n_part = 2,
    prop = 0.5
  )

  part_1$part
  part_1$best_part_info
  part_1$grid

  # Lets explore Grid object and presences and absences points
  plot(part_1$grid, col = gray.colors(20))
  points(part_1$part[c("x", "y")],
    col = rainbow(8)[part_1$part$.part],
    cex = 0.9,
    pch = c(1, 19)[part_1$part$pr_ab + 1]
  )
  expect_equal(class(part_1), "list")
  expect_equal(length(part_1), 3)
})

test_that("part_sband lat", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Let's practice with four latitudinal partition and only with presences
  single_spp <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
  part_2 <- part_sband(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    type = "lat",
    min_bands = 8,
    max_bands = 10,
    n_part = 8,
    prop = 0.5
  )

  part_2$part
  part_2$best_part_info
  part_2$grid
  expect_equal(class(part_2), "list")
  expect_equal(length(part_2), 3)
})


test_that("part_sband conditions for not finding a solution", {
  require(terra)
  require(dplyr)

  # Load datasets
  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Let's practice with four latitudinal partition and only with presences
  single_spp <- spp %>% dplyr::filter(species == "sp3", pr_ab == 1)

  expect_message(part_sband(
    env_layer = somevar,
    data = single_spp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    type = "lat",
    min_bands = 8,
    max_bands = 10,
    n_part = 8,
    prop = 0.5
  ))
})
