## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("class and lenght of gbm_t object", {
  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "cwd", "awc", "depth"
      ),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    )

  expect_equal(class(gbm_t), "list")
  expect_equal(length(gbm_t), 6)
})

test_that("test of 0-1 response argument", {
  data(abies)

  # We will partition the data with the k-fold method
  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  expect_error(gbm_t <-
    tune_gbm(
      data = abies2,
      response = "aet",
      predictors = c(
        "aet", "cwd", "awc", "depth"
      ),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    ))
})

test_that("test NULL in predictors_f", {
  data(abies)

  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
        "ppt_jja", "pH", "awc", "depth"
      ),
      predictors_f = NULL,
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    )

  expect_equal(class(gbm_t), "list")
  expect_equal(length(gbm_t), 6)
})

test_that("test if remove NAs rows works", {
  data(abies)

  # We will partition the data with the k-fold method
  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  # Insert NAs in rows 3 and 4 for response column.
  abies2[3:4, 1] <- NA

  expect_message(gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
        "ppt_jja", "pH", "awc", "depth"
      ),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    ))

  # Compare if the 2 NAs were removed
  testthat:::compare.numeric(nrow(abies2), nrow(gbm_t$data_ens))
})

test_that("test fit_formula", {
  data(abies)

  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  expect_message(gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "ppt_jja", "depth"
      ),
      predictors_f = c("landform"),
      fit_formula = formula("pr_ab ~ aet + ppt_jja + depth + landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    ))
})

test_that("grid = NULL ", {
  data(abies)

  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  expect_message(gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "awc", "depth"
      ),
      predictors_f = c("landform"),
      partition = ".part",
      grid = NULL,
      thr = "max_sens_spec",
      metric = "TSS"
    ))
})


test_that("missuse of grid ", {
  data(abies)

  abies2 <- part_random(
    data = abies %>% dplyr::group_by(pr_ab) %>%
      dplyr::slice_sample(prop = .2),
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 2)
  )

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50),
      # shrinkage = c(0.1, 0.5),
      n.minobsinnode = c(1, 3)
    )

  expect_error(gbm_t <-
    tune_gbm(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet", "awc", "depth"
      ),
      predictors_f = c("landform"),
      partition = ".part",
      grid = grid,
      thr = "max_sens_spec",
      metric = "TSS"
    ))
})
