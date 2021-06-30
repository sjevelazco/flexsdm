#' Merge model performance tables
#'
#' @param models
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #' data(abies)
#' abies
#'
#' # We will partition the data with the k-fold method
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
#'   thr = "max_sens_spec"
#' )
#'
#' # Build a generalized linear model using fit_glm
#'
#' glm_t1 <- fit_glm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
#'   poly = 0,
#'   inter_order = 0
#' )
#'
#' # Build a tuned random forest model using tune_raf
#'
#' tune_grid <-
#'   expand.grid(mtry = seq(1, 7, 1))
#'
#' rf_t <-
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
#'     thr = "max_sens_spec",
#'     metric = "TSS",
#'   )
#'}
sdm_summarize <- function(models) {

  if (data.class(models ) != "list") {
    stop("models must be a list object")
  }

  perf <-
    vector(mode = "list", length = length(models)) # initiate a list for the performance tables

  for (i in 1:length(models)) {
    perf[[i]] <- models[[i]]$performance
  }

  perf_tib <- rbindlist(perf, fill = TRUE) %>% as_tibble()

  return(perf_tib)
}
