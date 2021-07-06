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
  expect_true(all(class(t1)%in%c("tbl_df", "tbl", "data.frame")))

  # Function with background
  t1 <- sdm_eval(p = p, a = a, bg = backg)
  expect_true(all(class(t1)%in%c("tbl_df", "tbl", "data.frame")))

  # Function with >1000 presences and absences
  t1 <- sdm_eval(p = rep(p, 100), a = rep(a, 100), bg = backg)
  expect_true(all(class(t1)%in%c("tbl_df", "tbl", "data.frame")))

  # Test sensitivity threshold
  t1 <- sdm_eval(p = p, a = a, bg = backg, thr=c('sensitivity', sens = 0.5))
  expect_true(all(class(t1) %in% c("tbl_df", "tbl", "data.frame")))

  # test an error based on the misuse of threshold argument
  expect_error(sdm_eval(p = p, a = a, bg = backg, thr = "asdf"))

  # test an error based on the misuse of threshold argument
  expect_error(sdm_eval(p =p, a = NULL, thr ="max_fpb"))
})

