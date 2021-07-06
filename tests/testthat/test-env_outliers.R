test_that("test with sp > and < 15 occurrences", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp_ <- spp %>% dplyr::filter(species == "sp2")
  spp_ <- spp_ %>% mutate(idd = 1:nrow(spp_))

  outs_1 <- env_outliers(
    data = spp_,
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar
  )

  expect_equal(ncol(outs_1 %>% select(starts_with("."))), 7)
  expect_equal(sum(outs_1$.out_sum, na.rm = TRUE), 12)
})


test_that("test with dataset  with < 15 occurrences and only presence", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp_ <- spp %>% dplyr::filter(species == "sp2", pr_ab==1)
  spp_ <- spp_ %>% mutate(idd = 1:nrow(spp_))

  set.seed(1)
  spp_ <- spp_[sample(1:nrow(spp_), 10), ]
  outs_2 <- env_outliers(
    data = spp_,
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar
  )
  expect_equal(ncol(outs_2 %>% select(starts_with("."))), 7)
  expect_equal(sum(outs_2$.out_lof, na.rm = TRUE), 1)
})


test_that("test NA filtering ", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp
  spp_ <- spp %>% dplyr::filter(species == "sp2", pr_ab==1)
  spp_ <- spp_ %>% mutate(idd = 1:nrow(spp_))
  spp_$x[1] <- (-300000)
  outs_1 <- env_outliers(
    data = spp_,
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar
  )

  expect_equal(sum(colSums(is.na(outs_1 %>% select(starts_with("."))))), 7)
})


test_that("test with occurrences fewer than 6", {
  require(dplyr)
  require(terra)

  # Environmental variables
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)

  # Species occurrences
  data("spp")
  spp_ <- spp %>% dplyr::filter(species == "sp3", pr_ab==1)
  spp_ <- spp_ %>% mutate(idd = 1:nrow(spp_))
  spp_ <- spp_[1:5,]
  outs_1 <- env_outliers(
    data = spp_,
    pr_ab = "pr_ab",
    x = "x",
    y = "y",
    id = "idd",
    env_layer = somevar
  )

  outs_1.1 <- outs_1 %>% select(starts_with(".")) %>% na.omit()
  expect_true(all(colSums(outs_1.1==0)==4) )
})
