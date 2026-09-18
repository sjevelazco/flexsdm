skip_if_not_installed("terra")

test_that("errors when 'to' is not specified", {
  d <- data.frame(x = -74.1, y = 4.6)
  expect_error(occ_reproj(d, to = NULL), "A projection must be specified")
  # 'to' is also NULL by default, so calling without it should error too
  expect_error(occ_reproj(d), "A projection must be specified")
})

test_that("appends x_new/y_new columns by default (replace_cols = FALSE)", {
  d <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    x = c(-74.1, -73.9, -74.0),
    y = c(4.60, 4.70, 4.65)
  )

  out <- occ_reproj(d, to = "EPSG:3857")

  expect_s3_class(out, "tbl_df")
  expect_true(all(c("x", "y", "x_new", "y_new") %in% names(out)))
  expect_equal(nrow(out), nrow(d))
  # original columns are untouched
  expect_equal(out$x, d$x)
  expect_equal(out$y, d$y)
  # non-coordinate columns are preserved
  expect_equal(out$species, d$species)
})

test_that("reprojected coordinates differ from the originals for a real transform", {
  d <- data.frame(x = -74.1, y = 4.6)
  out <- occ_reproj(d, from = "EPSG:4326", to = "EPSG:3857")

  expect_true(is.numeric(out$x_new))
  expect_true(is.numeric(out$y_new))
  expect_false(isTRUE(all.equal(out$x_new, d$x)))
  expect_false(isTRUE(all.equal(out$y_new, d$y)))
})

test_that("replace_cols = TRUE overwrites the original x/y columns in place", {
  d <- data.frame(
    species = c("sp1", "sp2"),
    x = c(-74.1, -73.9),
    y = c(4.60, 4.70)
  )

  out <- occ_reproj(d, to = "EPSG:3857", replace_cols = TRUE)

  expect_equal(nrow(out), nrow(d))
  expect_named(out, names(d))
  expect_false(isTRUE(all.equal(out$x, d$x)))
  expect_false(isTRUE(all.equal(out$y, d$y)))
  expect_equal(out$species, d$species)
})

test_that("round-tripping the projection returns coordinates close to the original", {
  d <- data.frame(x = -74.1, y = 4.6)

  forward <- occ_reproj(d, x = "x", y = "y", from = "EPSG:4326", to = "EPSG:3857")
  back <- occ_reproj(
    data.frame(x = forward$x_new, y = forward$y_new),
    x = "x", y = "y", from = "EPSG:3857", to = "EPSG:4326"
  )

  expect_equal(back$x_new, d$x, tolerance = 1e-6)
  expect_equal(back$y_new, d$y, tolerance = 1e-6)
})

test_that("custom x/y column names are respected", {
  d <- data.frame(lon = c(-74.1, -73.9), lat = c(4.60, 4.70))

  out <- occ_reproj(d, x = "lon", y = "lat", to = "EPSG:3857")

  expect_true(all(c("lon", "lat", "x_new", "y_new") %in% names(out)))
  expect_equal(nrow(out), nrow(d))
})

test_that("works with a single row of data", {
  d <- data.frame(x = -74.1, y = 4.6)
  out <- occ_reproj(d, to = "EPSG:3857")

  expect_equal(nrow(out), 1)
  expect_true(all(c("x_new", "y_new") %in% names(out)))
})

test_that("works with tibble input, not just data.frame", {
  skip_if_not_installed("tibble")
  d <- tibble::tibble(x = c(-74.1, -73.9), y = c(4.60, 4.70))

  out <- occ_reproj(d, to = "EPSG:3857")

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), nrow(d))
})

test_that("default 'from' CRS is EPSG:4326 (WGS84)", {
  d <- data.frame(x = -74.1, y = 4.6)

  out_default <- occ_reproj(d, to = "EPSG:3857")
  out_explicit <- occ_reproj(d, from = "EPSG:4326", to = "EPSG:3857")

  expect_equal(out_default$x_new, out_explicit$x_new)
  expect_equal(out_default$y_new, out_explicit$y_new)
})
