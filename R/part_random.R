#' Conventional data partitioning methods
#'
#' @description This function provides different conventional (randomized, non-spatial) partitioning
#' methods based on cross validation folds (kfold, rep_kfold, and loocv), as well as bootstrap (boot)
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records
#' for a given species
#' @param pr_ab character. Column name of "data" with presences, presence-absence, or
#' pseudo-absence. Presences must be represented by 1 and absences by 0
#' @param method character. Vector with data partitioning method to be used:
#' \itemize{
#'   \item kfold: Random partitioning into k-folds for cross-validation. Usage
#'   part=c(method= 'kfold', folds='5'). 'folds' refers to the number of folds for data
#'   partitioning, it assumes value >=1. Usage method = c(method = "kfold", folds = 10).
#'   \item rep_kfold: Random partitioning into repeated k-folds for  cross-validation.
#'   Usage method = c(method = "rep_kfold", folds = 10, replicates=10). 'folds' refers to the
#'   number of folds for data partitioning, it assumes value >=1. 'replicate' refers to the
#'   number of replicates, it assumes a value >=1.
#'   \item loocv: Leave-one-out cross-validation (a.k.a. Jackknife). It is a special case of k-fold
#'   cross validation where the number of partitions is equal to the number of records.
#'   Usage method = c(method = "loocv").
#'   \item boot: Random bootstrap partitioning. Usage part=c(method='boot', replicates='2',
#'   proportion='0.7'). 'replicate' refers to the number of replicates, it assumes a value >=1.
#'   'proportion' refers to the proportion of occurrences used for model fitting, it assumes a
#'   value >0 and <=1. In this example proportion='0.7' mean that 70\% of data will be used for
#'   model training, while 30\% for model testing.
#'   }
#'
#' @return
#' A tibble object with information used in the 'data' argument and additional columns named .part
#' containing the partition groups. The rep_kfold and boot method will return as many ".part"
#' columns as replicated defined. For the rest of the methods, a single .part column is returned.
#' For kfold, rep_kfold, and loocv methods partition, groups are defined by integers.
#' On the contrary, for boot method, the partition groups are defined by the characters
#' 'train' and 'test'.
#'
#' @references
#' \itemize{
#' \item Fielding, A. H., & Bell, J. F. (1997). A review of methods for the assessment of
#' prediction errors in conservation presence/absence models. Environmental Conservation,
#' 24(1), 38-49. https://doi.org/10.1017/S0376892997000088
#' }
#' #' @export
#'
#' @importFrom dplyr %>% group_by mutate n summarise filter select slice_sample full_join left_join tibble
#'
#' @seealso \code{\link{part_sblock}}, \code{\link{part_senv}}, \code{\link{sample_pseudoabs}}, \code{\link{sample_background}}
#'
#' @examples
#' \dontrun{
#' data("abies")
#' abies$partition <- NULL
#' abies <- tibble(abies)
#'
#' # K-fold method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 10)
#' )
#' abies2
#'
#' # Repeated K-fold method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 10, replicates = 10)
#' )
#' abies2
#'
#' # Leave-one-out cross-validation (loocv) method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "loocv")
#' )
#' abies2
#'
#' # Bootstrap method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "boot", replicates = 50, proportion = 0.7)
#' )
#' abies2
#' abies2$.part1 %>% table() # Note that for this method .partX columns have train and test words.
#' }
#'
part_random <- function(data, pr_ab, method = NULL) {
  if (!method[1] %in% c("kfold", "rep_kfold", "loocv", "boot")) {
    stop("method argument was missused, available methods area 'kfold', 'rep-kfold', 'loocv', or 'boot'")
  }
  .part <- BOOT1 <- BOOT2 <- boot <- NULL

  # kfold
  if (method["method"] == "kfold") {
    data <- data %>%
      dplyr::group_by(!!as.symbol(pr_ab)) %>%
      dplyr::mutate(.part = sample(rep(1:method["folds"], length.out = dplyr::n()))) %>%
      dplyr::group_by()
  }

  # rep_kfold
  if (method["method"] == "rep_kfold") {
    for (i in 1:method["replicates"]) {
      cname <- paste0(".part", i)
      data <-
        data %>%
        dplyr::group_by(!!as.symbol(pr_ab)) %>%
        dplyr::mutate(rep_kfold = sample(rep(1:method["folds"],
          length.out = dplyr::n()
        ))) %>%
        dplyr::group_by()
      colnames(data)[colnames(data) == "rep_kfold"] <- cname
    }
  }

  # loocv
  if (method["method"] == "loocv") {
    data <- data %>%
      dplyr::group_by(!!as.symbol(pr_ab)) %>%
      dplyr::mutate(.part = 1:dplyr::n()) %>%
      dplyr::group_by()
    filt <- data %>%
      dplyr::group_by(!!as.symbol(pr_ab)) %>%
      dplyr::summarise(max = max(.part))
    filtmi <- filt %>% dplyr::filter(max == min(max))
    if (nrow(filtmi) > 1) {
      filtmi <- filtmi[1, ]
    }
    filt <- data$.part > filtmi$max
    if (sum(filt) > 0) {
      data$.part[filt] <- rep(1:filtmi$max, length.out = sum(filt))
    }
  }

  # BOOOTSRAP
  if (method["method"] == "boot") {
    reps <- as.numeric(method["replicates"])
    prop <- as.numeric(method["proportion"])
    prop2 <- 1 - prop
    data <- data %>% dplyr::group_by()
    data <- data %>% dplyr::mutate(IDBOOT = 1:nrow(data))
    for (i in 1:reps) {
      data2 <- data %>% dplyr::select({{ pr_ab }}, "IDBOOT")
      data_train <- data2 %>%
        dplyr::group_by(!!as.symbol(pr_ab)) %>%
        dplyr::slice_sample(prop = prop) %>%
        dplyr::mutate(BOOT1 = "train")
      data_ttest <- data2 %>%
        dplyr::group_by(!!as.symbol(pr_ab)) %>%
        dplyr::slice_sample(prop = prop2) %>%
        dplyr::mutate(BOOT2 = "test")

      data2 <- dplyr::full_join(data_train, data_ttest, by = c("pr_ab", "IDBOOT")) %>%
        dplyr::mutate(boot = paste(BOOT1, BOOT2, sep = "-")) %>%
        dplyr::select(-BOOT1, -BOOT2) %>%
        dplyr::mutate(boot = gsub(paste(c("-NA", "NA-"), collapse = "|"), "", boot))

      data <- dplyr::left_join(data, data2, by = c("pr_ab", "IDBOOT"))
      colnames(data)[colnames(data) == "boot"] <- paste0(".part", i)
    }
    data$IDBOOT <- NULL
  }

  return(dplyr::tibble(data))
}
