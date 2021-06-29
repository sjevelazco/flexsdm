test_that("multiplication works", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  glm_t1 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    poly = 0,
    inter_order = 0
  )

  expect_equal(class(glm_t1), "list")

  # testing with polynomial and interaction term
  glm_t2 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    poly = 2,
    inter_order = 1
  )

  expect_equal(class(glm_t2), "list")

  # Using repeated k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 3, replicates = 3)
  )

  glm_t3 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
    poly = 3,
    inter_order = 2
  )
  expect_equal(class(glm_t3), "list")
})
