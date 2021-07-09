#' Create directories for saving the outputs of the flexsdm workflow
#'
#' @param main_dir character. Directory path containing main folder for saving model inputs and outputs
#' @param proj_dir logical. If TRUE, function will create a folder for predictors for different regions or time periods used to project models. (default FALSE).
#' @param save_results logical. If TRUE, function will create a folder for model results. (default FALSE).
#' @param models vector. Vector of model types that will be used. (default NULL)
#' @param ensemble vector. Vector of methods used to ensemble different models. (default NULL)
#' @param proj_folders vector. Vector of folder names for future scenarios/different regions/time periods to save model projections models.
#' @param msdm logical. If TRUE, function will create a folder for the output of two methods for correcting SDM overprediction (see msdm_priori and msdm_posteriori). (default FALSE)
#'
#' @return None
#'
#'
#' @export
#'
#' @examples
#'
#' # Create directory for SDM workflow
#' my_dir <- file.path(getwd(), 'flexsdm_example')
#' dir.create(my_dir)
#'
#' # Implement sdm_directory
#' sdm_directory(main_dir = my_dir,
#' proj_dir = TRUE,
#' save_results = TRUE,
#' models = c('GLM', 'GAM', 'RANDOM_FOREST'),
#' ensemble = c('MEAN', "MEANSUP"),
#' msdm = TRUE)
#'
#'
#'
#'
sdm_directory <-
  function(main_dir,
           proj_dir = FALSE,
           save_results = FALSE,
           models = NULL,
           ensemble = NULL,
           proj_folders = NULL,
           msdm = FALSE) {
    if (!dir.exists(main_dir)) {
      stop("main directory does not exist")
    }

    dir.create(file.path(main_dir, 'Occurrences'))
    dir.create(file.path(main_dir, 'Predictors'))

    if (proj_dir != FALSE) {
      dir.create(file.path(main_dir, 'Projection'))
      dir.create(file.path(main_dir, 'Projection/Projection_PCA'))
    }

    if (save_results != FALSE) {
      dir.create(file.path(main_dir, 'Results'))

      if (!is.null(models)) {
        dir.create(file.path(main_dir, 'Results/Algorithm'))
        model_folders <-
          file.path(main_dir, 'Results/Algorithm', models)
        sapply(model_folders, dir.create)

        if (proj_dir != FALSE) {
          current_dir <- file.path(model_folders, 'Current')
          sapply(current_dir, dir.create)

          future_dir <- file.path(model_folders, 'Future')
          sapply(future_dir, dir.create)

          if (!is.null(proj_folders)) {
            for (i in 1:length(future_dir)) {
              future_dir2 <- file.path(future_dir[[i]], proj_folders)
              sapply(future_dir2, dir.create)
            }

          }
        }
      }

      if (!is.null(ensemble)) {
        dir.create(file.path(main_dir, 'Results/Ensemble'))
        ensemble_folders <-
          file.path(main_dir, 'Results/Ensemble', ensemble)
        sapply(ensemble_folders, dir.create)

        if (proj_dir != FALSE) {
          current_dir <- file.path(ensemble_folders, 'Current')
          sapply(current_dir, dir.create)

          future_dir <- file.path(ensemble_folders, 'Future')
          sapply(future_dir, dir.create)

        }

        if (!is.null(proj_folders)) {
          for (i in 1:length(future_dir)) {
            future_dir2 <- file.path(future_dir[[i]], proj_folders)
            sapply(future_dir2, dir.create)
          }
        }
      }

      if (msdm != FALSE) {
        dir.create(file.path(main_dir, 'Results/M_SDM'))
        dir.create(file.path(main_dir, 'Results/M_SDM/msdm_priori'))
        dir.create(file.path(main_dir, 'Results/M_SDM/msdm_posteriori'))
      }
    }




  }
