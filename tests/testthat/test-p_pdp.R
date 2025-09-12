library(terra)
library(dplyr)

somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar) # environmental data
names(somevar) <- c("aet", "cwd", "tmx", "tmn")
cat <- system.file("external/clusters.shp", package = "flexsdm")
cat <- terra::vect(cat)
cat$clusters <- paste0("c", cat$clusters)
cat <- terra::rasterize(cat, somevar, field = "clusters")
somevar <- c(somevar, cat)


data(abies)

set.seed(123)
abies2 <- abies %>%
  select(x, y, pr_ab) %>%
  group_by(pr_ab) %>%
  dplyr::slice_sample(n = 70)

set.seed(210)
abies2 <- sdm_extract(abies2,
  x = "x",
  y = "y",
  env_layer = somevar
) %>%
  part_random(
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )


test_that("test p_pdp with continuous and factor and gam", {
  m_ <- fit_gam(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"),
    predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec"), k = 3
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and gau", {
  m_ <- fit_gau(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and glm", {
  m_ <- fit_glm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and gbm", {
  m_ <- fit_gbm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and max", {
  m_ <- fit_max(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and net", {
  m_ <- fit_net(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and raf", {
  m_ <- fit_raf(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})

test_that("test p_pdp with continuous and factor and svm", {
  m_ <- fit_svm(
    data = abies2, response = "pr_ab", predictors = c("aet", "cwd", "tmx", "tmn"), predictors_f = "clusters", partition = ".part", thr = c("max_sens_spec")
  )

  asdf <- p_pdp(model = m_$model, training_data = abies2)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 4)

  asdf <- p_pdp(model = m_$model, training_data = abies2, predictors = c("aet", "cwd"))
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(length(asdf$patches$plots), 1)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  asdf <- p_pdp(model = m_$model, training_data = abies2, resid = TRUE, resolution = 5)
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 5)

  # Partial depence plot for training and projection condition found in a projection area
  expect_equal(asdf$theme$line$colour, "black")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar, resolution = 50,
    colorl = c("#CC00FF", "#CCFF00")
  )
  expect_equal(class(asdf)[1], "patchwork")
  expect_equal(nrow(asdf$data), 52)

  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray"
  )
  expect_equal(class(asdf)[1], "patchwork")
  asdf <- p_pdp(
    model = m_$model, training_data = abies2, projection_data = somevar,
    colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
    theme = ggplot2::theme_dark()
  )
  asdf$data
  expect_equal(class(asdf)[1], "patchwork")
})
