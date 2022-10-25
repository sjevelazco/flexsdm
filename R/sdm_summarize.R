#' Merge model performance tables
#'
#' @param models list of one or more models fitted with fit_ or tune_ functions, or a fit_ensemble output, a esm_ family function output. A list a single or several models fitted with some of fit_ or tune_ functions or object returned by the \code{\link{fit_ensemble}} function. Usage models = list(mod1, mod2, mod3)
#'
#' @return
#'
#' Combined model performance table for all input models. Models fit with tune will include model performance for the best hyperparameters.
#'
#' @export
#'
#' @importFrom dplyr bind_rows relocate tibble
#'
#' @examples
#' \dontrun{
#' data(abies)
#' abies
#'
#' # In this example we will partition the data using the k-fold method
#'
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' # Build a generalized additive model using fit_gam
#'
#' gam_t1 <- fit_gam(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
#' )
#'
#' gam_t1$performance
#'
#' # Build a generalized linear model using fit_glm
#'
#' glm_t1 <- fit_glm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   poly = 0,
#'   inter_order = 0
#' )
#'
#' glm_t1$performance
#'
#' # Build a tuned random forest model using tune_raf
#'
#' tune_grid <-
#'   expand.grid(mtry = seq(1, 7, 1))
#'
#' rf_t1 <-
#'   tune_raf(
#'     data = abies2,
#'     response = "pr_ab",
#'     predictors = c(
#'       "aet", "cwd", "tmin", "ppt_djf",
#'       "ppt_jja", "pH", "awc", "depth"
#'     ),
#'     predictors_f = c("landform"),
#'     partition = ".part",
#'     grid = tune_grid,
#'     thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'     metric = "TSS",
#'   )
#'
#' rf_t1$performance
#'
#' # Merge sdm performance tables
#'
#' merge_df <- sdm_summarize(models = list(gam_t1, glm_t1, rf_t1))
#'
#' merge_df
#' }
sdm_summarize <- function(models) {
  . <- model_ID <- model <- IMAE_sd <- NULL
  if (data.class(models) != "list") {
    stop("models must be a list object")
  }

  # if more than 1 model is provided
  if (length(models) > 1) {
    # list of each model's performance table
    perf <- lapply(models, function(x) {
      x$performance
    })

    # add unique idea for each model
    perf <-
      Map(cbind, perf, model_ID = 1:length(perf))

    # bind rows and move model_ID to first column
    perf_tib <-
      dplyr::bind_rows(perf) %>%
      dplyr::relocate(model_ID, .before = model) %>%
      dplyr::tibble()
  } else {
    perf_tib <- models[[1]]$performance
    perf_tib$model_ID <- 1
  }

  perf_tib <- perf_tib %>%
    dplyr::relocate( names(
      dplyr::select(perf_tib, model_ID:IMAE_sd)
    ))
  return(perf_tib)
}
