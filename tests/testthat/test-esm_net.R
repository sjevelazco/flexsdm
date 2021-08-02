test_that("ESM Artificial Neural Networks", {
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
  esm_net_t1 <- esm_net(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmin", "ppt_djf"),
    partition = ".part",
    thr = NULL
  )

  expect_equal(length(esm_net_t1), 3)
  expect_equal(class(esm_net_t1), "list")
})
