test_that("multiplication works", {
  data("abies_db")

  # Using k-fold partition method
  abies_db2 <- part_random(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  svm_t1 <- fit_svm(
    data = abies_db2,
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
  abies_db2 <- part_random(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 3, proportion = 0.7)
  )

  svm_t2 <- fit_svm(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    fit_formula = NULL
  )

  expect_equal(class(svm_t2), "list")
})
