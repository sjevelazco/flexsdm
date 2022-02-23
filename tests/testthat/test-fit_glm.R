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
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    poly = 0,
    inter_order = 0
  )

  expect_equal(class(glm_t1), "list")

  # testing with polynomial and interaction term
  glm_t2 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
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
    predictors = c("ppt_jja", "pH"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    poly = 3,
    inter_order = 2
  )
  expect_equal(class(glm_t3), "list")


  # Does the function work without predictors_f?
  glm_t3 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_equal(class(glm_t3), "list")

  # What about no predictors? Does not work
  expect_error(fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  ))
})


test_that("test glm with NA, no factor variable and using formula", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  abies2$aet[2:10] <- NA
  glm_t1 <- fit_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = as.formula("pr_ab ~ aet + ppt_jja + pH + awc + depth")
  )

  expect_equal(class(glm_t1), "list")
})


test_that("test select_var argument", {
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
    select_pred = TRUE,
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    poly = 2,
    inter_order = 1
  )

  expect_equal(class(glm_t1), "list")
})
