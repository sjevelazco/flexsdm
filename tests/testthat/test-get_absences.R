test_that("errors when required columns are missing", {
  d <- data.frame(species = c("sp1", "sp2"), x = 1:2, y = 3:4)

  expect_error(get_absences(d, x = "lon"), "must contain the columns")
  expect_error(get_absences(d, y = "lat"), "must contain the columns")
  expect_error(get_absences(d, species = "sp"), "must contain the columns")
})

test_that("returns a tibble with the expected column names", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("species", "x", "y", "pr_ab"))
})

test_that("default pr_ab_name is 'pr_ab' and values are coded 1/0", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d)

  expect_true(all(out$pr_ab %in% c(0, 1)))

  # within the sp1 block, presences (1) should be exactly the rows whose
  # coordinates match the original sp1 occurrences
  block_sp1 <- out[out$species == "sp1", ]
  sp1_coords <- d[d$species == "sp1", c("x", "y")]
  expect_equal(
    sort(block_sp1$pr_ab[
      paste(block_sp1$x, block_sp1$y) %in% paste(sp1_coords$x, sp1_coords$y)
    ]),
    rep(1, nrow(sp1_coords))
  )
})

test_that("custom pr_ab_name is respected (regression test for the dead-argument bug)", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d, pr_ab_name = "occ")

  expect_true("occ" %in% names(out))
  expect_false("pr_ab" %in% names(out))
  expect_false("prba" %in% names(out))
  expect_true(all(out$occ %in% c(0, 1)))
})

test_that("with target_species = NULL, builds one block per unique species", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d)

  n_species <- length(unique(d$species))
  expect_equal(nrow(out), nrow(d) * n_species)
  expect_setequal(unique(out$species), unique(d$species))

  # each species block should have exactly one presence row per original
  # occurrence of that species, and absences for every other row
  for (s in unique(d$species)) {
    block <- out[out$species == s, ]
    expect_equal(nrow(block), nrow(d))
    expect_equal(sum(block$pr_ab == 1), sum(d$species == s))
    expect_equal(sum(block$pr_ab == 0), sum(d$species != s))
  }
})

test_that("target_species restricts output to the requested species only", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d, target_species = "sp1")

  expect_equal(nrow(out), nrow(d))
  expect_equal(unique(out$species), "sp1")
  expect_equal(sum(out$pr_ab == 1), sum(d$species == "sp1"))
  expect_equal(sum(out$pr_ab == 0), sum(d$species != "sp1"))
})

test_that("target_species with multiple species stacks one block per target", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
    x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
    y = c(4.60, 4.65, 4.70, 4.72, 4.68)
  )

  out <- get_absences(d, target_species = c("sp1", "sp2"))

  expect_equal(nrow(out), nrow(d) * 2)
  expect_setequal(unique(out$species), c("sp1", "sp2"))
})

test_that("a target species not present in the data still produces an all-absence block", {
  d <- data.frame(
    species = c("sp1", "sp1", "sp2", "sp2"),
    x = c(-74.1, -74.2, -73.9, -73.8),
    y = c(4.60, 4.65, 4.70, 4.72)
  )

  out <- get_absences(d, target_species = "sp_unseen")

  expect_equal(nrow(out), nrow(d))
  expect_true(all(out$species == "sp_unseen"))
  expect_true(all(out$pr_ab == 0))
})

test_that("only species, x, y (and pr_ab) columns are retained in output", {
  d <- data.frame(
    species = c("sp1", "sp2"),
    x = c(-74.1, -73.9),
    y = c(4.60, 4.70),
    elevation = c(120, 340),
    notes = c("a", "b")
  )

  out <- get_absences(d)

  expect_setequal(names(out), c("species", "x", "y", "pr_ab"))
  expect_false("elevation" %in% names(out))
  expect_false("notes" %in% names(out))
})

test_that("custom x/y/species column names are respected", {
  d <- data.frame(
    sp_name = c("sp1", "sp1", "sp2"),
    lon = c(-74.1, -74.2, -73.9),
    lat = c(4.60, 4.65, 4.70)
  )

  out <- get_absences(d, x = "lon", y = "lat", species = "sp_name")

  expect_named(out, c("sp_name", "lon", "lat", "pr_ab"))
  expect_equal(nrow(out), nrow(d) * length(unique(d$sp_name)))
})

test_that("works with tibble input, not just data.frame", {
  skip_if_not_installed("tibble")
  d <- tibble::tibble(
    species = c("sp1", "sp1", "sp2"),
    x = c(-74.1, -74.2, -73.9),
    y = c(4.60, 4.65, 4.70)
  )

  out <- get_absences(d)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), nrow(d) * length(unique(d$species)))
})

test_that("unique species are sorted when target_species = NULL", {
  d <- data.frame(
    species = c("sp3", "sp1", "sp2"),
    x = c(-74.0, -74.1, -73.9),
    y = c(4.68, 4.60, 4.70)
  )

  out <- get_absences(d)

  # blocks are stacked in the order of `sp`, so the first appearance of
  # each species in the output should follow sorted order
  expect_equal(unique(out$species), sort(unique(d$species)))
})
