test_that("part_senv", {
  require(terra)

  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Select a species
  spp1 <- spp %>% dplyr::filter(species == "sp1")

  # Test with presences and absences
  part1 <- part_senv(
    env_layer = somevar,
    data = spp1,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_n_groups = 2,
    max_n_groups = 10,
    prop = 0.2
  )

  expect_equal(class(part1), "list")

  # Only with presences
  spp1 <- spp1 %>% dplyr::filter(pr_ab == 1)
  part2 <- part_senv(
    env_layer = somevar,
    data = spp1,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_n_groups = 2,
    max_n_groups = 20,
    prop = 0.2
  )

  expect_equal(class(part2), "list")
})


test_that("misuse of arguments", {
  skip_on_cran()
  require(terra)

  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  # Select a species
  spp1 <- spp %>% dplyr::filter(species == "sp1")
  spp1$pr_ab[1:2] <- 3
  # Test with presences and absences
  expect_error(part_senv(
    env_layer = somevar,
    data = spp1,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_n_groups = 2,
    max_n_groups = 10,
    prop = 0.2
  ))
})

test_that("include_coords = FALSE is independent of user column names", {
  skip_on_cran()
  require(terra)

  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)

  spp1 <- spp %>% dplyr::filter(species == "sp1")

  # Columns are renamed to literal "pr_ab"/"x"/"y" on entry, so the
  # include_coords = FALSE exclusion must use those literal names. Excluding
  # by the user's argument values kept pr_ab and the coordinates in the
  # k-means input whenever the input columns were not already named
  # pr_ab/x/y, changing the partition.
  set.seed(1)
  part_default <- part_senv(
    env_layer = somevar,
    data = spp1,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    min_n_groups = 2,
    max_n_groups = 10,
    prop = 0.2,
    include_coords = FALSE
  )

  spp2 <- spp1 %>% dplyr::rename(occ = pr_ab, lon = x, lat = y)
  set.seed(1)
  part_renamed <- part_senv(
    env_layer = somevar,
    data = spp2,
    x = "lon",
    y = "lat",
    pr_ab = "occ",
    min_n_groups = 2,
    max_n_groups = 10,
    prop = 0.2,
    include_coords = FALSE
  )

  expect_identical(part_default$part$.part, part_renamed$part$.part)
})
