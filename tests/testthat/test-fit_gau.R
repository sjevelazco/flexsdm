test_that("multiplication works", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies2

  bg <- abies2
  bg$pr_ab <- 0


  gaup_t1 <- fit_gau(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    background = bg,
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_equal(class(gaup_t1), "list")

  # Using bootstrap partition method and only with presence-absence
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 5, proportion = 0.7)
  )

  gaup_t2 <- fit_gau(
    data = abies2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c(type = c("lpt", "max_sens_spec", "sensitivity"), sens = "0.8")
  )

  expect_equal(class(gaup_t2), "list")

  # Does the function work without predictors_f?
  gaup_t3 <- fit_gau(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_equal(class(gaup_t3), "list")

  # What about no predictors? Does not work
  expect_error(fit_gau(
    data = abies2,
    response = "pr_ab",
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  ))
})


test_that("test gau with NA, no factor variable and using formula", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  abies2$aet[2:10] <- NA
  gaup_t1 <- fit_gau(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_equal(class(gaup_t1), "list")
})


