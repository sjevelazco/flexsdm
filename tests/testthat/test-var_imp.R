test_that("var_imp", {
  require(dplyr)

  data(abies)
  data(backg)
  set.seed(0)
  abies <- abies %>%
    dplyr::group_by(pr_ab) %>%
    dplyr::slice_sample(prop = .2)
  set.seed(0)
  backg <- backg %>%
    dplyr::group_by(pr_ab) %>%
    dplyr::slice_sample(prop = .2)

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  max_t1 <- fit_max(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin"),
    predictors_f = c("landform"),
    partition = ".part",
    background = backg2,
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    clamp = TRUE,
    classes = "default",
    pred_type = "cloglog",
    regmult = 1
  )


  net_t1 <- fit_net(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  svm_f1 <- fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_message(
    v_ip <-
      var_imp(
        data = abies2,
        response = "pr_ab",
        predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin", "landform"),
        models = list(max_t1, net_t1, svm_f1),
        clamp = TRUE,
        pred_type = "cloglog",
        thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
        n_sim = 50,
        n_cores = 5
      )
  )
  expect_length(v_ip, 13)
})
