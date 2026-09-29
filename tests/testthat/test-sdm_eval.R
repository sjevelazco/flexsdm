test_that("fist test", {
  require(dplyr)

  set.seed(0)
  p <- rnorm(50, mean = 0.7, sd = 0.3) %>% abs()
  p[p > 1] <- 1
  p[p < 0] <- 0

  set.seed(0)
  a <- rnorm(50, mean = 0.3, sd = 0.2) %>% abs()
  a[a > 1] <- 1
  a[a < 0] <- 0

  set.seed(0)
  backg <- rnorm(1000, mean = 0.4, sd = 0.4) %>% abs()
  backg[backg > 1] <- 1
  backg[backg < 0] <- 0

  # Function use without threshold specification
  t1 <- sdm_eval(p = p, a = a)
  expect_true(all(class(t1) %in% c("tbl_df", "tbl", "data.frame")))

  # Function with background
  t1 <- sdm_eval(p = p, a = a, bg = backg)
  expect_true(all(class(t1) %in% c("tbl_df", "tbl", "data.frame")))

  # Function with >1000 presences and absences
  t1 <- sdm_eval(p = rep(p, 100), a = rep(a, 100), bg = backg)
  expect_true(all(class(t1) %in% c("tbl_df", "tbl", "data.frame")))

  # Test sensitivity threshold
  t1 <- sdm_eval(p = p, a = a, bg = backg, thr = c("sensitivity", sens = 0.5))
  expect_true(all(class(t1) %in% c("tbl_df", "tbl", "data.frame")))

  # test an error based on the misuse of threshold argument
  expect_error(sdm_eval(p = p, a = a, bg = backg, thr = "asdf"))

  # test an error based on the misuse of threshold argument
  expect_error(sdm_eval(p = p, a = NULL, thr = "max_fpb"))
})

test_that("KAPPA is Cohen's kappa", {
  # Worked confusion matrix: the max_sens_spec threshold here is 0.9,
  # which gives tp=40, fn=10, fp=5, tn=45
  p <- c(rep(0.9, 40), rep(0.1, 10))
  a <- c(rep(0.9, 5), rep(0.1, 45))

  e <- sdm_eval(p = p, a = a, thr = "max_sens_spec")

  # By hand: Pr(a) = (40+45)/100 = 0.85
  # Pr(e) = ((40+5)*(40+10) + (10+45)*(5+45)) / 100^2 = (2250 + 2750)/10000 = 0.5
  # kappa = (0.85 - 0.5)/(1 - 0.5) = 0.7
  expect_equal(e$KAPPA, 0.7)

  # An uninformative model (identical scores for presences and absences) must
  # score 0, and a perfect model must score 1
  e0 <- sdm_eval(p = rep(0.8, 50), a = rep(0.8, 50), thr = "max_sens_spec")
  expect_equal(e0$KAPPA, 0)
  e1 <- sdm_eval(p = rep(0.9, 50), a = rep(0.1, 50), thr = "max_sens_spec")
  expect_equal(e1$KAPPA, 1)
})
