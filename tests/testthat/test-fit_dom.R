test_that("fit_dom works with k-fold partition and continuous + categorical predictors", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  dom_t1 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
  )

  expect_equal(class(dom_t1), "list")
  expect_named(
    dom_t1,
    c("model", "predictors", "performance", "performance_part", "data_ens")
  )

  # model$domain must contain only presence rows for the modeled predictors
  expect_true(all(c("aet", "ppt_jja", "pH", "awc", "depth", "landform") %in%
    names(dom_t1$model$domain)))

  # one performance row per threshold type requested
  expect_equal(nrow(dom_t1$performance), 3)
  expect_true(all(c("max_sens_spec", "equal_sens_spec", "max_sorensen") %in%
    dom_t1$performance$threshold))
  expect_true(all(dom_t1$performance$TPR_mean >= 0 & dom_t1$performance$TPR_mean <= 1))
  expect_true(all(dom_t1$performance$AUC_mean >= 0 & dom_t1$performance$AUC_mean <= 1))

  # performance_part has one row per replica x partition x threshold
  expect_equal(nrow(dom_t1$performance_part), 3 * 3)
  expect_true(all(c("replica", "partition", "model", "threshold") %in%
    names(dom_t1$performance_part)))

  # data_ens carries one predicted row per test observation, usable by fit_ensemble
  expect_true(all(c("rnames", "replicates", "part", "pr_ab", "pred") %in%
    names(dom_t1$data_ens)))
})


test_that("fit_dom works without predictors_f (continuous predictors only)", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  dom_t2 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  expect_equal(class(dom_t2), "list")
  expect_false("landform" %in% names(dom_t2$model$domain))
  expect_equal(nrow(dom_t2$performance), 1)
})


test_that("fit_dom works with repeated k-fold partitioning", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 3, replicates = 3)
  )

  dom_t3 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("max_sens_spec", "equal_sens_spec")
  )

  expect_equal(class(dom_t3), "list")
  # 3 replicates x 3 folds x 2 thresholds
  expect_equal(nrow(dom_t3$performance_part), 3 * 3 * 2)
  expect_equal(length(unique(dom_t3$performance_part$replica)), 3)
})


test_that("fit_dom works with multiple threshold types, including sensitivity", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  dom_t4 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = c("lpt", "max_jaccard", "sensitivity", sens = "0.8")
  )

  expect_equal(class(dom_t4), "list")
  expect_true(all(c("lpt", "max_jaccard", "sensitivity") %in%
    dom_t4$performance$threshold))
})


test_that("fit_dom uses all default thresholds when thr = NULL", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  dom_t5 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = ".part",
    thr = NULL
  )

  expect_equal(class(dom_t5), "list")
  expect_gt(nrow(dom_t5$performance), 1)
})


test_that("fit_dom with partition = NULL returns only the stored presence model", {
  data("abies")

  dom_t6 <- fit_dom(
    data = abies,
    response = "pr_ab",
    predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
    predictors_f = c("landform"),
    partition = NULL,
    thr = c("max_sens_spec")
  )

  expect_equal(class(dom_t6), "list")
  expect_named(dom_t6, "model")
  expect_true(all(abies$pr_ab[rownames(dom_t6$model$domain) %>% as.numeric()] == 1))
})


test_that("fit_dom removes rows with NAs in predictors and reports it", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )
  abies2$aet[1:10] <- NA

  expect_message(
    dom_t7 <- fit_dom(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
      predictors_f = c("landform"),
      partition = ".part",
      thr = c("max_sens_spec")
    ),
    "rows were excluded"
  )

  expect_equal(class(dom_t7), "list")
})


test_that("fit_dom errors when no predictors are provided", {
  data("abies")

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  expect_error(
    fit_dom(
      data = abies2,
      response = "pr_ab",
      predictors_f = c("landform"),
      partition = ".part",
      thr = c("max_sens_spec")
    )
  )
})


test_that("fit_dom model object can be used by sdm_predict", {
  skip_if_not_installed("terra")
  data("abies")
  somevar <- terra::rast(system.file("external/somevar.tif", package = "flexsdm"))
  names(somevar) <- c("aet", "cwd", "tmx", "tmn")

  abies2 <- abies %>%
    dplyr::select(x, y, pr_ab) %>%
    dplyr::group_by(pr_ab) %>%
    dplyr::slice_sample(prop = 0.5) %>%
    dplyr::ungroup()
  abies2 <- sdm_extract(abies2, x = "x", y = "y", env_layer = somevar)
  abies2 <- part_random(
    data = abies2,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  dom_t8 <- fit_dom(
    data = abies2,
    response = "pr_ab",
    predictors = c("aet", "cwd", "tmx", "tmn"),
    predictors_f = NULL,
    partition = ".part",
    thr = c("max_sens_spec")
  )

  pred <- sdm_predict(
    models = dom_t8,
    pred = somevar,
    thr = "max_sens_spec",
    con_thr = FALSE,
    predict_area = NULL
  )

  expect_true(is.list(pred))
  expect_true(methods::is(pred$dom$dom, "SpatRaster"))
})
