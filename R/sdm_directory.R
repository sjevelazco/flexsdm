#' Create directories for saving the outputs of the flexsdm
#'
#' @description This function assists in creating a directory system with different sub-folders to
#' assist in the organisation of the modelling process outputs.
#'
#' @param main_dir character. Directory path containing main folder for saving model inputs and
#' outputs. There is no default: the path must be provided by the user (e.g., in examples use
#' a sub-folder of \code{tempdir()}). The folder is created if it does not exist.
#' @param projections vector. Vector of folder names for future scenarios/different regions/time
#' periods to save model projections output.
#' @param algorithm vector. Vector of model names that will be used. Usage algorithm =
#' c(gam, tune_max, tune_net, esm_glm). If  "all" is used the function creates folders with all
#' algorithms available in flexsdm . i.e. 'gam', 'gau', 'gbm', 'glm', 'max', 'net', 'raf', and
#' 'svm'. Default NULL
#' @param calibration_area logical. If TRUE, function creates a folder in 1_Inputs for storing
#' calibration area. Default TRUE
#' @param ensemble vector. Vector of methods used to ensemble different models.
#' Usage ensemble = c("mean", "meanthr"). Default NULL
#' @param threshold logical. If TRUE sub-folders "/1_con", "/2_bin" will be created within each
#' algorithm and/or ensemble folder. Used for storing continuous and binarized models separately.
#' Default FALSE
#' @param return_vector logical. If TRUE function returns a vector with the path to all folders.
#' Default TRUE
#'
#' @return A character vector with the paths to created folders
#'
#' @details
#' The sdm_directory function assists in saving workflow outputs by
#' creating folders (directories) based on user specifications, such as choice
#' of algorithms, ensemble methods, and model projections to new geographic
#' regions or periods. The function first creates two folders within the
#' user-specified project folder, one for model inputs (1_Inputs) and one for
#' model outputs (2_Outputs). Within 1_Inputs, there are three sub-folders for
#' users to store model inputs: 1_Occurrences, 2_Predictors, and 3_Calibration_area.
#' If a user chooses to include projections in their modeling framework, the
#' 2_Projections subfolder is created within the 2_Predictors folder to store
#' the environmental data for the projection scenarios provided to the
#' "projections" argument. Additionally, sdm_directory offers users enhanced
#' flexibility in saving their modeling outputs, giving them offers users the
#' option to save results from any modeling and ensemble technique presented in flexsdm
#'
#' @export
#'
#' @importFrom dplyr %>% tibble
#'
#' @examples
#' # Folders are created in a temporary directory; in a real analysis use
#' # the path of your own project folder
#' my_project <- file.path(tempdir(), "my_project_name")
#' # my_project <- file.path(getwd(), "flexsdm_results") # create folders in the working directory
#' dirs_1 <- sdm_directory(
#'   main_dir = my_project,
#'   projections = NULL,
#'   calibration_area = TRUE,
#'   algorithm = c("gam", "tune_max"),
#'   ensemble = c("mean", "meanthr"),
#'   threshold = FALSE,
#'   return_vector = TRUE
#' )
#' dirs_1
#' unlink(my_project, recursive = TRUE) # remove the folders created in this example
#'
#' # With folders for model projections and for thresholded outputs
#' dirs_2 <- sdm_directory(
#'   main_dir = my_project,
#'   projections = c("cnrm_rpc8.5_2050", "cnrm_rpc4.5_2050"),
#'   calibration_area = TRUE,
#'   algorithm = "all",
#'   ensemble = c("mean", "meanthr"),
#'   threshold = TRUE
#' )
#' head(dirs_2)
#' unlink(my_project, recursive = TRUE)
sdm_directory <-
  function(
    main_dir,
    projections = NULL,
    calibration_area = TRUE,
    algorithm = NULL,
    ensemble = NULL,
    threshold = FALSE,
    return_vector = TRUE
  ) {
    . <- NULL
    if (missing(main_dir) || is.null(main_dir)) {
      stop(
        "'main_dir' must be provided: path of the folder where the project directories will be created"
      )
    }
    if (!dir.exists(main_dir)) {
      dir.create(main_dir, recursive = TRUE)
    }

    # Create input directories
    dir.create(
      file.path(main_dir, "1_Inputs/1_Occurrences"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(main_dir, "1_Inputs/2_Predictors/1_Current"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (isTRUE(calibration_area)) {
      dir.create(
        file.path(main_dir, "1_Inputs/3_Calibration_area"),
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    # Create output directories
    out_dir <- file.path(main_dir, "2_Outputs")
    dir.create(
      file.path(out_dir, "0_Model_performance"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(out_dir, "1_Current"),
      recursive = TRUE,
      showWarnings = FALSE
    )

    # Handle projections
    out_dir_proj <- NULL
    if (!is.null(projections)) {
      # Input Projections
      in_proj_dir <- file.path(main_dir, "1_Inputs/2_Predictors/2_Projection")
      dir.create(in_proj_dir, recursive = TRUE, showWarnings = FALSE)
      sapply(
        file.path(in_proj_dir, projections),
        dir.create,
        showWarnings = FALSE
      )

      # Output Projections
      out_proj_dir <- file.path(out_dir, "2_Projection")
      dir.create(out_proj_dir, recursive = TRUE, showWarnings = FALSE)
      out_dir_proj <- file.path(out_proj_dir, projections)
      sapply(out_dir_proj, dir.create, showWarnings = FALSE)
    }

    # Process algorithm 'all' keyword
    if (!is.null(algorithm) && any("all" == algorithm)) {
      algorithm <- c("gam", "gau", "gbm", "glm", "max", "net", "raf", "svm")
    }

    # Helper function to create Algorithm/Ensemble structures efficiently
    create_structure <- function(
      parent_dirs,
      type_name,
      items,
      threshold_flag
    ) {
      if (is.null(items)) {
        return()
      }
      for (pd in parent_dirs) {
        base_path <- file.path(pd, type_name)
        dir.create(base_path, showWarnings = FALSE)

        item_paths <- file.path(base_path, items)
        sapply(item_paths, dir.create, showWarnings = FALSE)

        if (threshold_flag) {
          sapply(
            file.path(item_paths, "1_con"),
            dir.create,
            showWarnings = FALSE
          )
          sapply(
            file.path(item_paths, "2_bin"),
            dir.create,
            showWarnings = FALSE
          )
        }
      }
    }

    # Identify all parent directories where model outputs should be stored
    parents <- file.path(out_dir, "1_Current")
    if (!is.null(out_dir_proj)) {
      parents <- c(parents, out_dir_proj)
    }

    # Create Algorithm and Ensemble directories across all contexts
    create_structure(parents, "Algorithm", algorithm, threshold)
    create_structure(parents, "Ensemble", ensemble, threshold)

    message("Directories were created in:\n", main_dir)

    if (return_vector) {
      results <- list.dirs(main_dir)
      # results <- gsub(getwd(), '.', results)
      return(results)
    }
  }
