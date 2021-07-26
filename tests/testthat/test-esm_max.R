test_that("ESM maximum entropy", {
  data("abies")
  data("backg")
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
    method = c(method = "rep_kfold", folds = 5, replicates = 5)
  )
  abies2

  set.seed(10)
  backg2 <- backg %>%
    na.omit() %>%
    group_by(pr_ab) %>%
    dplyr::slice_sample(n = 100) %>%
    group_by()

  backg2 <- part_random(
    data = backg2,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 5, replicates = 5)
  )
  backg2

  # Without threshold specification and with kfold
  esm_max_t1 <- esm_max(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmin", "ppt_djf"),
    partition = ".part",
    thr = NULL,
    background = backg2,
    clamp = TRUE,
    classes = "default",
    pred_type = "cloglog",
    regmult = 1
  )

  expect_equal(length(esm_max_t1), 3)
  expect_equal(class(esm_max_t1), "list")
})
