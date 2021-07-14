test_that("mask method", {
  require(terra)
  require(dplyr)
  require(rgdal)
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
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})

test_that("mcp method", {
  require(terra)
  require(dplyr)
  require(rgdal)
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
    method = "mcp",
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})

test_that("bmcp method", {
  require(terra)
  require(dplyr)
  require(rgdal)
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

test_that("bmcp method", {
  require(terra)
  require(dplyr)
  require(rgdal)
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
  require(rgdal)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)

  single_spp <- single_spp %>% mutate(groups = ifelse(x > 150000, "a", "b"))

  # buffer method
  ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("bmcp", width = 40000), groups = "groups"
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})


test_that("mask method", {
  require(terra)
  require(dplyr)
  require(rgdal)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  # clusters <- terra::vect(clusters)
  # To see if it will be convert to SpatVector object
  clusters <- rgdal::readOGR(clusters)

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
    method = c("mask", clusters, "clusters"),
  )

  expect_equal(class(ca_1)[1], "SpatVector")
})


## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("missuse method argument", {
  require(terra)
  require(dplyr)
  require(rgdal)
  data("spp")
  clusters <- system.file("external/clusters.shp", package = "flexsdm")
  clusters <- terra::vect(clusters)

  single_spp <-
    spp %>%
    dplyr::filter(species == "sp1") %>%
    dplyr::filter(pr_ab == 1) %>%
    dplyr::select(-pr_ab)


  # buffer method
  expect_error(ca_1 <- calib_area(
    data = single_spp,
    x = "x",
    y = "y",
    method = c("maskclustersbu1"),
  ))
})
