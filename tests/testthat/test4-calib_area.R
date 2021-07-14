test_that("multiplication works", {
  require(terra)
  require(dplyr)
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
