test_that("multiplication works", {
  require(dplyr)
  require(terra)
  require(ggplot2)

  # Envirnomental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp
  spp1 <- spp %>% dplyr::filter(species == "sp3")
  spp1 <- spp1 %>% mutate(idd = 1:nrow(spp1))

  outs_1 <- env_outliers(
    data = spp1,
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar
  )
  expect_equal(ncol(outs_1 %>% select(starts_with('.'))), 7)
  expect_equal(sum(outs_1$.out_sum, na.rm = TRUE), 7)
})

