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
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  gbm_t <-
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
    )

  expect_equal(class(gbm_t), "list")
  expect_equal(length(gbm_t), 5)

  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test of 0-1 response argument", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  expect_error(gbm_t <-
                 tune_gbm(
                   data = abies2,
                   response = "aet",
                   predictors = c(
                     "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
                     "ppt_jja", "pH", "awc", "depth"
                   ),
                   predictors_f = c("landform"),
                   partition = ".part",
                   grid = tune_grid,
                   thr = "max_sens_spec",
                   metric = "TSS"
                 )
               )

  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test of non character predictor", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  expect_error(gbm_t <-
                 tune_gbm(
                   data = abies2,
                   response = "pr_ab",
                   predictors = c(
                     "landform"
                   ),
                   predictors_f = c("landform"),
                   partition = ".part",
                   grid = tune_grid,
                   thr = "max_sens_spec",
                   metric = "TSS"
                 )
  )

  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test NULL in predictors_f", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
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
  expect_equal(length(gbm_t), 5)

  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test if remove NAs rows works", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  # Insert NAs in rows 3 and 4 for response column.
  abies2[3:4,1] <- NA

  expect_error(gbm_t <-
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
    )
  )

  # Compare if the 2 NAs were removed
  testthat:::compare.numeric(nrow(abies2), nrow(gbm_t$data_ens))

  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test fit_formula", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  expect_message(gbm_t <-
                 tune_gbm(
                   data = abies2,
                   response = "pr_ab",
                   predictors = c(
                     "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
                     "ppt_jja", "pH", "awc", "depth"
                   ),
                   predictors_f = c("landform"),
                   fit_formula = pr_ab ~ aet + ppt_jja + depth + landform,
                   partition = ".part",
                   grid = tune_grid,
                   thr = "max_sens_spec",
                   metric = "TSS"
                 )
  )


  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test grid is NULL", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

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
                     grid = NULL,
                     thr = "max_sens_spec",
                     metric = "TSS"
                   )
  )


  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)


test_that("test fit_formula", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  expect_message(gbm_t <-
                   tune_gbm(
                     data = abies2,
                     response = "pr_ab",
                     predictors = c(
                       "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
                       "ppt_jja", "pH", "awc", "depth"
                     ),
                     predictors_f = c("landform"),
                     fit_formula = pr_ab ~ aet + ppt_jja + depth + landform,
                     partition = ".part",
                     grid = tune_grid,
                     thr = "max_sens_spec",
                     metric = "TSS"
                   )
  )


  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)

test_that("test grid == character", {

  data(abies)
  abies

  # We will partition the data with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  tune_grid <-
    expand.grid(
      n.trees = c(20, 50, 100),
      shrinkage = c(0.1, 0.5, 1),
      n.minobsinnode = c(1, 3, 5, 7, 9)
    )

  # Test if grid == character
  tune_grid2 <- as.character(tune_grid)

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
                     grid = tune_grid2,
                     thr = "max_sens_spec",
                     metric = "TSS"
                   )
  )


  # Outputs
  gbm_t$model
  gbm_t$predictors
  gbm_t$performance
  gbm_t$hyper_performance
  gbm_t$data_ens

  # Graphical exploration of performance of each hyper-parameter setting
  require(ggplot2)
  pg <- position_dodge(width = 0.5)
  ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
                                      TSS_mean,
                                      col = factor(shrinkage)
  )) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
                  width = 0.2, position = pg
    ) +
    geom_point(position = pg) +
    geom_line(
      data = gbm_t$tune_performance,
      aes(as.numeric(factor(n.minobsinnode)),
          TSS_mean,
          col = factor(shrinkage)
      ), position = pg
    ) +
    facet_wrap(. ~ n.trees) +
    theme(legend.position = "bottom")
}
)
