## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("test example tune_max", {
  require(maxnet)
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

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(
      regmult = seq(0.1, 3, 0.5),
      classes = c("l", "lq")
    )

  expect_message(max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = c("landform"),
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    ))
  expect_length(max_t, 5)
})


test_that("test NULL predictors_f and NULL grid", {
  require(maxnet)
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

  # We will partition the data and background with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  abies2

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  backg2

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(
      regmult = seq(0.1, 3, 1),
      classes = c("l", "lq")
    )

  max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = NULL,
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    )
  expect_length(max_t, 5)

  max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = NULL,
      background = backg2,
      partition = ".part",
      grid = NULL,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    )
  expect_length(max_t, 5)
})

test_that("test data with NA and without background", {
  require(maxnet)
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

  abies[c(1, 50, 100), 2:ncol(abies)] <- NA
  backg[c(1, 50, 100), 2:ncol(backg)] <- NA

  # We will partition the data and background with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  abies2

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  backg2

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(
      regmult = seq(0.1, 3, 1),
      classes = c("l", "lq")
    )

  max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = NULL,
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    )
  expect_length(max_t, 5)

  max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = NULL,
      background = NULL,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    )
  expect_length(max_t, 5)

})

test_that("test fit only with presences and background", {
  require(maxnet)
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

  # We will partition the data and background with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  # Only presences
  abies2 <- abies2 %>% dplyr::filter(pr_ab==1)

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  backg2

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(
      regmult = seq(0.1, 3, 1),
      classes = c("l", "lq")
    )

  max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = "landform",
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = FALSE,
      pred_type = "cloglog"
    )
  expect_length(max_t, 5)

})

test_that("test background argument names not match
          lack of hyperparameter", {
  require(maxnet)
  require(dplyr)

  data(abies)
  abies

  data(backg)
  backg

  # We will partition the data and background with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  abies2

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )
  backg2

  # Change 1st column name to test error message
  backg2 <- backg2 %>%
    dplyr::rename(pr_ab_RNMD = pr_ab)

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(
      regmult = seq(0.1, 31),
      classes = c("l", "lq")
    )

  expect_error(max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = c("landform"),
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    ))

  gridtest <-
    expand.grid(
      regmult = seq(0.1, 31),
      classes = c("l", "lq")
    )

  expect_error(max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = c("landform"),
      background = backg2,
      partition = ".part",
      grid = gridtest[-1],
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    ))

  # Groups doesn't match
  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 3)
  )

  expect_error(max_t <-
    tune_max(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "tmin"),
      predictors_f = c("landform"),
      background = backg2,
      partition = ".part",
      grid = gridtest,
      thr = "max_sens_spec",
      metric = "TSS",
      clamp = TRUE,
      pred_type = "cloglog"
    ))
})
