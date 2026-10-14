test_that("sdm_directory", {
  # All folders are created inside tempdir() (never in the working directory)
  main_dir_1 <- file.path(tempdir(), "flexsdm_results_test")

  dirs_1 <- sdm_directory(
    main_dir = main_dir_1,
    projections = NULL,
    calibration_area = TRUE,
    algorithm = c("gam", "tune_max"),
    ensemble = c("mean", "meanthr"),
    threshold = FALSE,
    return_vector = TRUE
  )
  expect_equal(class(dirs_1), "character")
  expect_equal(length(dirs_1), 15)
  expect_true(all(dir.exists(dirs_1)))
  unlink(main_dir_1, recursive = TRUE)

  main_dir_2 <- file.path(tempdir(), "my_project_name")
  dirs_2 <- sdm_directory(
    main_dir = main_dir_2,
    projections = c(
      "cnrm_rpc8.5_2050",
      "cnrm_rpc4.5_2050"
    ),
    calibration_area = TRUE,
    algorithm = "all",
    ensemble = c("mean", "meanthr"),
    threshold = TRUE
  )
  expect_equal(class(dirs_2), "character")
  expect_equal(length(dirs_2), 111)
  unlink(main_dir_2, recursive = TRUE)
})


test_that("sdm_directory requires an explicit main_dir", {
  # the function must never fall back to writing in the working directory
  expect_error(sdm_directory(), "main_dir")
  expect_error(sdm_directory(main_dir = NULL), "main_dir")
})
