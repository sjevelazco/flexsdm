#' Merge model performance tables
#'
#' @param models
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' models <- list(m_1, m_2, m_3)
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
