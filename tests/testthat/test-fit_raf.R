test_that("test raf", {
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

  # Does the function work without predictors_f?
  raf_t3 <- fit_raf(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(raf_t3), "list")

  # What about no predictors? Does not work
  expect_error(fit_raf(
    data = abies2,
    response = "pr_ab",
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = NULL
  ))
})


test_that("test raf with NA, no factor variable and using formula", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  abies2$aet[2:10] <- NA
  rf_t1 <- fit_raf(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = as.formula("pr_ab ~ aet + ppt_jja + pH + awc + depth")
  )

  expect_equal(class(rf_t1), "list")
})
