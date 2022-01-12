test_that("plot res", {
  f <- system.file("external/somevar.tif", package = "flexsdm")
  r <- terra::rast(f)
  r <- r$CFP_1
  expect_equal(class(plot_res(r, res_mult = 100)), "list")
})
