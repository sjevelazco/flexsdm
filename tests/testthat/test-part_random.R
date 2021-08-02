test_that("part_random", {
  data("abies")
  abies$partition <- NULL
  abies <- tibble(abies)

  # K-fold method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 10)
  )
  expect_true(all(unique(abies2$.part) %in% 1:10))

  # Repeated K-fold method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 10, replicates = 10)
  )
  abies2
  expect_true(all(unique(abies2$.part1) %in% 1:10))
  expect_equal(ncol(abies2 %>% dplyr::select(starts_with(".part"))), 10)

  # Leave-one-out cross-validation (loocv) method
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "loocv")
  )
  expect_true(length(unique(abies2 %>% dplyr::filter(pr_ab == 0) %>% pull(.part))) == 251)

  # Bootstrap method
  set.seed(10)
  abies2 <- part_random(
    data = abies,
    pr_ab = "pr_ab",
    method = c(method = "boot", replicates = 50, proportion = 0.7)
  )
  expect_equal(length(unique(abies2$.part1)), 4)
  expect_equal((abies2 %>% dplyr::select(dplyr::starts_with(".part"))) %>% ncol(), 50)
})
