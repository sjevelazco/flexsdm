test_that("sdm_directory", {
  require(dplyr)
  # require(sf)

  # Implement sdm_directory without specific path and project name

  dirs_1 <- sdm_directory(
    main_dir = NULL,
    projections = NULL,
    calibration_area = TRUE,
    algorithm = c("gam", "tune_max"),
    ensemble = c("mean", "meanthr"),
    threshold = FALSE,
    return_vector = TRUE
  )
  expect_equal(class(dirs_1), "character")
  expect_equal(length(dirs_1), 15)
  unlink(file.path(getwd(), "flexsdm_results"), recursive = TRUE)

  # Implement sdm_directory without specific path and project name
  dirs_2 <- sdm_directory(
    main_dir = tempdir() %>% file.path(., "my_project_name"),
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
  unlink(tempdir() %>% file.path(., "my_project_name"), recursive = TRUE)
})
