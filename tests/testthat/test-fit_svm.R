test_that("multiplication works", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  svm_t1 <- fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(svm_t1), "list")

  # Using bootstrap partition method and only with presence-absence
  # and get performance for several method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 3, proportion = 0.7)
  )

  svm_t2 <- fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(svm_t2), "list")

  # Does the function work without predictors_f?
  svm_t3 <- fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(svm_t3), "list")

  # What about no predictors? Does not work
  expect_error(fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    fit_formula = NULL
  ))

})
