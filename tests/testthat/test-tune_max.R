## %######################################################%##
#                                                          #
####          Set of tests for testing errors           ####
#                                                          #
## %######################################################%##

test_that("test example tune_max", {
  require(maxnet)

  data(abies)
  abies

  data(backg)
  backg

  # We will partition the data and background with the k-fold method

  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )
  abies2

  backg2 <- part_random(
    data = backg,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 5)
  )
  backg2

  # pr_ab columns is species presence and absences (i.e. the response variable)
  # from aet to landform are the predictors variables (landform is a qualitative variable)

  # Hyper-parameter values for tuning
  gridtest <-
    expand.grid(regmult = seq(0.1, 3, 0.5),
                classes = c("l", "lq", "lqh", "lqhp", "lqhpt"))

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
    )
  )
}
)
