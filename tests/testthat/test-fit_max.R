test_that("multiplication works", {
  data("abies_db")

  # Using k-fold partition method
  abies_db2 <- part_classical(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  # generating background data
  bg <- abies_db2
  bg$pr_ab <- 0

  max_t1 <- fit_max(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    background = bg,
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen")
  )

  expect_equal(class(max_t1), "list")

  # Using bootstrap partition method
  abies_db2 <- part_classical(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 5, proportion = 0.7)
  )

  # generating background data
  bg <- abies_db2
  bg$pr_ab <- 0

  max_t2 <- fit_max(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    background = bg,
    thr = c(type = c("lpt", "max_sens_spec", "sensitivity"), sens = "0.8")
  )

  expect_equal(class(max_t2), "list")
})
