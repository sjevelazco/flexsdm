test_that("ESM generalized linear model ", {
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
    method = c(method = "kfold", folds = 3)
  )
  abies2

  # Without threshold specification and with kfold
  esm_gam_t1 <- esm_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmin", "ppt_djf"),
    partition = ".part",
    thr = NULL
  )

  expect_equal(class(esm_gam_t1), "list")
  expect_equal(length(esm_gam_t1), 3)
})
