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
suppressMessages(
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
)

suppressMessages(
  svm_f1 <- fit_svm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )
)


test_that("sdm_varimp", {
  v_ip <-
    sdm_varimp(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "cwd", "tmin", "landform"),
      models = list(max_t1, svm_f1),
      clamp = TRUE,
      pred_type = "cloglog",
      thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
      n_sim = 10,
      n_cores = 1
    )
  expect_length(v_ip, 13)
})

test_that("esemble", {
  vmensemble <- flexsdm::fit_ensemble(
    models = list(svm_f1, max_t1),
    ens_method = "mean",
    thr = NULL,
    thr_model = "max_sens_spec",
    metric = "TSS"
  )

  v_ip <- sdm_varimp(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH"),
    models = vmensemble,
    clamp = TRUE,
    pred_type = "cloglog",
    thr = "max_sens_spec",
    n_sim = 10,
    n_cores = 1
  )

  expect_true(all(c("svm", "max") %in% v_ip$model))
})

test_that("sdm_varimp for esm", {
  data("abies")
  require(dplyr)

  # Using k-fold partition method
  set.seed(10)
  abies2 <- abies %>%
    na.omit() %>%
    group_by(pr_ab) %>%
    dplyr::slice_sample(n = 10) %>%
    group_by()

  abies2 <- part_random(
    data = abies2,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 3, replicates = 5)
  )
  abies2

  # Without threshold specification and with kfold
  esm_glm_t1 <- esm_glm(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmin", "ppt_djf"),
    partition = ".part",
    thr = NULL,
    poly = 1,
    inter_order = 1
  )


  v_ip <-
    sdm_varimp(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "ppt_jja", "pH"),
      models = esm_glm_t1,
      clamp = TRUE,
      pred_type = "cloglog",
      thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
      n_sim = 10,
      n_cores = 1
    )
  expect_true(all(grepl("esm", v_ip$model)))
})
