test_that("multiplication works", {
  data("abies_db")

  # Using k-fold partition method
  abies_db2 <- part_random(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies_db2

  gam_t1 <- fit_gam(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec"
  )

  expect_equal(class(gam_t1), "list")


  # Using our own formula
  require(gam)
  gam_t2 <- fit_gam(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec",
    fit_formula = stats::formula(pr_ab ~ s(aet, df = 4) +
      s(ppt_jja, df = 3) +
      s(pH, df = 3) + landform)
  )

  expect_equal(class(gam_t2), "list")

  # Using repeated k-fold partition method
  abies_db2 <- part_random(
    data = abies_db,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 5, replicates = 5)
  )

  gam_t3 <- fit_gam(
    data = abies_db2,
    response = "pr_ab",
    predictors = c("ppt_jja", "pH", "awc"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = "max_sens_spec"
  )

  expect_equal(class(gam_t3), "list")
})
