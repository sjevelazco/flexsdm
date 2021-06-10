test_that("multiplication works", {

  data("abies_db")

  # Using k-fold partition method
  abies_db2 <- part(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies_db2

  gbm_t1 <- fit_gbm(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen")
  )

  expect_equal(class(gbm_t1), "list")

  # Using bootstrap partition method
  abies_db2 <- part(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 10, proportion = 0.7)
  )
  abies_db2

  gbm_t2 <- fit_gbm(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec"
  )

  expect_equal(class(gbm_t2), "list")
})
