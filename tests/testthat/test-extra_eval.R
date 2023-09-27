require(dplyr)
require(dplyr)
require(terra)

data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)
names(somevar) <- c("aet", "cwd", "tmx", "tmn")


spp$species %>% unique()
sp <- spp %>%
  dplyr::filter(species == "sp3", pr_ab == 1) %>%
  dplyr::select(x, y, pr_ab)

# Calibration area based on some criterion such as dispersal ability
ca <- calib_area(sp,
  x = "x", y = "y",
  method = c("bmcp", width = 50000),
  crs = crs(somevar)
)

# plot(somevar[[1]])
# points(sp)
# plot(ca, add = T)

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

# Merge presences and absences databases to get a complete calibration data
sp_pa <- dplyr::bind_rows(sp, psa)
sp_pa

# Get environmental condition of calibration area
sp_pa_2 <- sdm_extract(data = sp_pa, x = "x", y = "y", env_layer = somevar)



test_that("extra_eval one and two cores", {
  # Measure degree of extrapolation based on Mahalanobis and
  # for a projection area based on a SpatRaster object
  extr <-
    extra_eval(
      training_data = sp_pa_2,
      projection_data = somevar,
      pr_ab = "pr_ab",
      n_cores = 1,
      aggreg_factor = 1,
      metric = "mahalanobis"
    )
  expect_equal(class(extr)[1], "SpatRaster")

  # Based on euclidean distance
  extr <-
    extra_eval(
      training_data = sp_pa_2,
      projection_data = somevar,
      pr_ab = "pr_ab",
      n_cores = 1,
      aggreg_factor = 1,
      metric = "euclidean"
    )
  expect_equal(class(extr)[1], "SpatRaster")
})

test_that("extra_eval with uni_comb argument", {
  # Measure degree of extrapolation based on Mahalanobis
  extr <-
    extra_eval(
      training_data = sp_pa_2,
      projection_data = somevar,
      pr_ab = "pr_ab",
      n_cores = 1,
      aggreg_factor = 1,
      metric = "mahalanobis",
      univar_comb = TRUE
    )
  expect_equal(class(extr)[1], "SpatRaster")
  expect_equal(names(extr), c("extrapolation", "uni_comb"))
})


test_that("extra_eval based on tibble object", {
  # Based on euclidean distance and dataframe
  extr <-
    extra_eval(
      training_data = sp_pa_2,
      projection_data = as.data.frame(somevar) %>% as_tibble(),
      pr_ab = "pr_ab",
      n_cores = 1,
      aggreg_factor = 1,
      metric = "euclidean"
    )
  expect_equal(class(extr)[1], "tbl_df")
})

test_that("extra_eval based on tibble object with uni_comb argument", {
  # Based on euclidean distance and dataframe
  extr <-
    extra_eval(
      training_data = sp_pa_2,
      projection_data = as.data.frame(somevar) %>% as_tibble(),
      pr_ab = "pr_ab",
      n_cores = 1,
      aggreg_factor = 1,
      metric = "euclidean",
      univar_comb = TRUE
    )
  expect_equal(class(extr)[1], "tbl_df")
  expect_equal(names(extr)[1:2],c("extrapolation", "univar_comb"))
})


test_that("extra_eval wrong use", {
  expect_error(extra_eval(
    training_data = sp_pa_2,
    projection_data = somevar,
    n_cores = 1,
    aggreg_factor = 3
  ))

  expect_error(extra_eval(
    training_data = sp_pa_2,
    projection_data = as.data.frame(somevar) %>% as_tibble(),
    pr_ab = "pr_ab",
    n_cores = 1,
    aggreg_factor = 1,
    metric = "euclidean12321"
  ))
})
