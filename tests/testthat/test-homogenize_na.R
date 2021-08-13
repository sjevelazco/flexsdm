test_that("test homogenize_na", {
  require(dplyr)
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)[[1]]
  v2 <- system.file("external/suit_time_step.tif", package = "flexsdm")
  v2 <- terra::rast(v2)[[1]]
  somevar <- terra::rast(list(v2, somevar))
  somevar2 <- homogenize_na(somevar)
  sum1 <- terra::global(is.na(somevar[[1]]), sum)
  sum2 <- terra::global(is.na(somevar2[[1]]), sum)
  expect_true(sum1 < sum2)
})
