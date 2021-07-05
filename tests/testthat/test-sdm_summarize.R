test_that("multiplication works", {

  data(abies)

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # Build a generalized additive model using fit_gam

  gam_t1 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  # Build a generalized linear model using fit_glm

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

  # Build a tuned random forest model using tune_raf

  tune_grid <-
    expand.grid(mtry = seq(1, 7, 1))

  rf_t1 <-
    tune_raf(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin", "ppt_djf",
                     "ppt_jja", "pH", "awc", "depth"),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
      metric = "TSS",
    )


  # Merge sdm performance tables

  merge_df <- sdm_summarize(models = list(gam_t1, glm_t1, rf_t1))

  expect_equal(data.class(merge_df), "data.frame")


  # only provide one model to sdm_summarize
  # must specify as list(model)

  merge_df2 <- sdm_summarize(models = list(gam_t1))

  expect_equal(data.class(merge_df2), "data.frame")

})
