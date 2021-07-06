test_that("multiplication works", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  rf_t1 <- fit_raf(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(rf_t1), "list")

  # Using bootstrap partition method and only with presence-absence
  # and get performance for several methods
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 3, proportion = 0.7)
  )

  rf_t2 <- fit_raf(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(rf_t2), "list")
})
