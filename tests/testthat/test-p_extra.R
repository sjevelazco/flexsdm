test_that("test p_extra", {
  require(dplyr)
  require(terra)
  require(ggplot2)

  data(spp)
  f <- system.file("external/somevar.tif", package = "flexsdm")
  somevar <- terra::rast(f)
  names(somevar) <- c("aet", "cwd", "tmx", "tmn")

  spp$species %>% unique()
  sp <- spp %>%
    dplyr::filter(species == "sp2", pr_ab == 1) %>%
    dplyr::select(x, y, pr_ab)

  # Calibration area based on some criterion such as dispersal ability
  ca <- calib_area(sp, x = "x", y = "y", method = c("buffer", width = 50000), crs = crs(somevar))

  # Sampling pseudo-absences
  set.seed(10)
  psa <- sample_pseudoabs(
    data = sp,
    x = "x",
    y = "y",
    n = nrow(sp) * 2, # selecting number of pseudo-absence points twice number of presences
    method = "random",
    rlayer = somevar,
    calibarea = ca
  )

  # Merge presences and abasences databases to get a complete calibration data
  sp_pa <- dplyr::bind_rows(sp, psa)
  sp_pa

  # Get environmental condition of calibration area
  sp_pa_2 <- sdm_extract(data = sp_pa, x = "x", y = "y", env_layer = somevar)
  sp_pa_2

  # Measure extrapolation based on calibration data (presence and pseudo-absences)
  # using SHAPE metric
  extr <-
    extra_eval(
      training_data = sp_pa_2, # change by training_data
      pr_ab = "pr_ab",
      projection_data = somevar, # change to projection_data
      n_cores = 1,
      aggreg_factor = 1
    )

  ##%######################################################%##
  ####            Explore extrapolation in the            ####
  ####        environmental and geographical space        ####
  ##%######################################################%##

  asfd <- p_extra(
    training_data = sp_pa_2,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    extra_suit_data = extr,
    projection_data = somevar,
    geo_space = TRUE,
    prop_points = 0.01
  )
  expect_equal(class(asfd), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asfd$data), 33207)


  asdf <- p_extra(
    training_data = sp_pa_2,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    extra_suit_data = extr,
    projection_data = somevar,
    predictors = c("tmn", "cwd"),
    geo_space = TRUE,
    prop_points = 0.05
  )
  expect_equal(class(asfd), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asfd$data), 33207)

  asdf <- p_extra(
    training_data = sp_pa_2,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    extra_suit_data = extr,
    projection_data = somevar,
    predictors = c("cwd", "tmx", "aet"),
    geo_space = TRUE,
    geo_position = "left",
    prop_points = 0.05,
    color_p = "white",
    alpha_p = 0.5,
    alpha_gradient = 0.2,
    color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
    theme = ggplot2::theme_dark()
  )

  expect_equal(class(asfd), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asfd$data), 33207)

  asdf <- p_extra(
    training_data = sp_pa_2,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    extra_suit_data = extr,
    projection_data = somevar,
    geo_space = TRUE,
    prop_points = 0.05,
    color_p = "white",
    alpha_p = 0.5,
    alpha_gradient = 0.2,
    color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
    theme = ggplot2::theme_dark()
  )

  expect_equal(class(asfd), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asfd$data), 33207)

  # Explore extrapolation only in the environmental space
  asdf <- p_extra(
    training_data = sp_pa_2,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    extra_suit_data = extr,
    projection_data = somevar,
    geo_space = FALSE,
    prop_points = 0.05,
    color_p = "black",
    color_gradient = c("#085CF8", "#65AF1E", "#F3CC1D", "#FC6A9B", "#D70500"),
    theme = ggplot2::theme_minimal()
  )

  expect_equal(class(asfd), c("patchwork", "gg", "ggplot"))
  expect_equal(nrow(asfd$data), 33207)
})
