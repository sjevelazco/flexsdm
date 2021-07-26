test_that("test homogenize_na", {
  somevar <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(somevar)[[c(1, 4)]]
  somevar2 <- homogenize_na(somevar)
  sum1 <- terra::global(is.na(somevar[[2]]), sum)
  sum2 <- terra::global(is.na(somevar2[[2]]), sum)
  expect_true(sum1 < sum2)
})
