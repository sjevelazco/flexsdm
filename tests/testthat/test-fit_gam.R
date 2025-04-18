test_that("multiplication works", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies2

  gam_t1 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec"
  )

  expect_equal(class(gam_t1), "list")


  # Using our own formula
  require(mgcv)
  gam_t2 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec",
    fit_formula = stats::formula(pr_ab ~ s(aet, k = 4) +
      s(ppt_jja, k = 3) +
      s(pH, k = 3) + landform)
  )

  expect_equal(class(gam_t2), "list")

  # Using repeated k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 5, replicates = 5)
  )

  gam_t3 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec"
  )

  expect_equal(class(gam_t3), "list")
})

test_that("test gam with NA, no factor variable and using formula", {
  data("abies")

  # Using k-fold partition method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  abies2$aet[2:10] <- NA
  gam_t1 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
    fit_formula = as.formula("pr_ab ~ aet + ppt_jja + pH + awc + depth")
  )

  expect_equal(class(gam_t1), "list")
})


test_that("fit_gam with error: number of data minor than parameters ", {
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
  expect_message(esm_gam_t1 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmin", "ppt_djf"),
    partition = ".part",
    thr = NULL,
    k = 10
  ))
  expect_equal(class(esm_gam_t1), "NULL")
})



test_that("test select_var argument", {
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies2

  gam_t3 <- fit_gam(
    data = abies2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    select_pred = TRUE,
    partition = ".part",
    thr = "max_sens_spec"
  )
  gam_t3

  expect_equal(class(gam_t3), "list")
})
