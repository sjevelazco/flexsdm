## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##
data(abies)

# We will partition the data with the k-fold method
abies2 <- part_random(
  data = abies %>% dplyr::group_by(pr_ab) %>%
    dplyr::slice_sample(prop = .2),
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 2)
)


test_that("tuen", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(200,400,600,800))

  raf_t <-
    tune_raf(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "awc", "depth"),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    )

  expect_equal(class(raf_t), "list")
  expect_equal(length(raf_t), 5)
})


test_that("tuen without ntree hyperparamenter", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1))

  raf_t <-
    tune_raf(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "cwd", "awc", "depth"),
      predictors_f = c("landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    )

  expect_equal(class(raf_t), "list")
  expect_equal(length(raf_t), 5)
})


test_that("test of 0-1 response argument", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(400,600))

  expect_error(
    raf_t <-
      tune_raf(
        data = abies2,
        response = "aet",
        predictors = c("aet", "cwd", "awc", "depth"),
        predictors_f = c("landform"),
        partition = ".part",
        grid = tune_grid,
        thr = "max_sens_spec",
        metric = "TSS"
      )
  )
})


test_that("test NULL in predictors_f", {
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(400,600))

  raf_t <-
    tune_raf(
      data = abies2,
      response = "pr_ab",
      predictors = c(
        "aet",
        "cwd",
        "tmin",
        "ppt_djf",
        "ppt_jja",
        "ppt_jja",
        "pH",
        "awc",
        "depth"
      ),
      predictors_f = NULL,
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS"
    )

  expect_equal(class(raf_t), "list")
  expect_equal(length(raf_t), 5)
})

test_that("test if remove NAs rows works", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(400,600))

  # Insert NAs in rows 3 and 4 for response column.
  abies2[3:4, 1] <- NA

  expect_message(
    raf_t <-
      tune_raf(
        data = abies2,
        response = "pr_ab",
        predictors = c(
          "aet",
          "cwd",
          "tmin",
          "ppt_djf",
          "ppt_jja",
          "ppt_jja",
          "pH",
          "awc",
          "depth"
        ),
        predictors_f = c("landform"),
        partition = ".part",
        grid = tune_grid,
        thr = "max_sens_spec",
        metric = "TSS"
      )
  )

  # Compare if the 2 NAs were removed
  testthat:::compare.numeric(nrow(abies2), nrow(raf_t$data_ens))
})

test_that("test fit_formula", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(400,600))

  raf_t <-
    tune_raf(
      data = abies2,
      response = "pr_ab",
      predictors = c("aet", "ppt_jja", "depth"),
      predictors_f = c("landform"),
      fit_formula = formula("pr_ab ~ aet + ppt_jja + depth + landform"),
      partition = ".part",
      grid = tune_grid,
      thr = "max_sens_spec",
      metric = "TSS",
      n_cores = 3
    )
  expect_equal(length(raf_t), 5)
})

test_that("grid = NULL ", {
  tune_grid <-
    expand.grid(mtry = seq(1, 4, 1),
                ntree  = c(400,600))
  expect_message(
    raf_t <-
      tune_raf(
        data = abies2,
        response = "pr_ab",
        predictors = c("aet", "awc", "depth"),
        predictors_f = c("landform"),
        partition = ".part",
        grid = NULL,
        thr = "max_sens_spec",
        metric = "TSS"
      )
  )
})


test_that("missuse of grid ", {
  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(mtrsdfy = seq(1, 4, 1),
                ntrdee  = c(400,600))

  expect_error(
    raf_t <-
      tune_raf(
        data = abies2,
        response = "pr_ab",
        predictors = c("aet", "awc", "depth"),
        predictors_f = c("landform"),
        partition = ".part",
        grid = tune_grid,
        thr = "max_sens_spec",
        metric = "TSS"
      )
  )
})
