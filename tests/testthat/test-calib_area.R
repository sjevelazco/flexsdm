test_that("buffer method", {
  require(terra)
  require(dplyr)
  require(rgeos)

  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # buffer method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("buffer", width = 40000),
    crs = NULL
  )

  expect_equal(class(ca_1)[1], "SpatVector")

  # # buffer method with crs
  # ca_1 <- calib_area(
  #   data = single_spp,
  #   x = "x",
  #   y = "y",
  #   method = c("buffer", width = 40000),
  #   crs = crs(clusters)
  # )
  #
  # expect_equal(class(ca_1)[1], "SpatVector")
})

test_that("mcp method", {
  require(terra)
  require(dplyr)
  require(rgeos)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # mcp method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = "mcp",
    crs = NULL
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})


test_that("mcp method with group", {
  require(terra)
  require(dplyr)
  require(rgeos)
  data("spp")

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)
  single_spp <- single_spp %>% mutate(groups = ifelse(x > 150000, "a", "b"))

  # mcp method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = "mcp",
    crs = NULL,
    groups = "groups"
  )

  expect_equal(class(ca_1)[1], "SpatVector")
  expect_equal(length(ca_1), 2)
})

test_that("bmcp method", {
  require(terra)
  require(dplyr)
  require(rgeos)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # buffer method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("bmcp", width = 40000), groups = NULL
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})

test_that("bmcp method with groups", {
  require(terra)
  require(dplyr)
  require(rgeos)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)
  single_spp <- single_spp %>% mutate(groups = ifelse(x > 150000, "a", "b"))


  # bmcp method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("bmcp", width = 40000),
    groups = "groups",
    crs = NULL
  )

  expect_equal(class(ca_1)[1], "SpatVector")
  expect_equal(length(ca_1), 2)
})


# test_that("mask method", {
#   require(terra)
#   require(dplyr)
#   require(rgeos)
#   data("spp")
#   clusters <- system.file("external/clusters.shp", package = "flexsdm")
#   clusters <- terra::vect(clusters)
#   # To see if it will be convert to SpatVector object
#   clusters <- as(clusters, "Spatial")
#
#   single_spp <-
#     spp %>%
#     dplyr::filter(species == "sp1") %>%
#     dplyr::filter(pr_ab == 1) %>%
#     dplyr::select(-pr_ab)
#
#
#   # buffer method
#   ca_1 <- calib_area(
#     data = single_spp,
#     x = "x",
#     y = "y",
#     method = c("mask", clusters, "clusters"),
#   )
#
#   expect_equal(class(ca_1)[1], "SpatVector")
# })


## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("missuse method argument", {
  require(terra)
  require(dplyr)
  require(rgeos)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # buffer method
  expect_error(calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("maskclustersbu1")
  ))

  expect_error(calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("bmcp")
  ))

  expect_error(calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("mask")
  ))
})
