#'
#' Data partitioning for training and testing models
#'
#' @description This function provides different conventional partition methods based in folds (kfold, rep_kfold, and loocv), and bootstrap (boot)
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records for a given species
#' @param pr_ab character. Column name of "data" with presences, presence-absence, or pseudo-absence. Presences must be represented by 1 and absences by 0
#' @param bg_data data.frame. Data frames with background points.
#' @param bg_a character. Column name of "bg_data" with absences. Background must be repented by 0
#' @param method character. Vector with partition method to be used:
#' \itemize{
#'   \item kfold: Random partition in k-fold cross-validation. Usage part=c(method= 'kfold', folds='5'). 'folds' refers to the number of folds for data partitioning, it assumes value >=1. Usage method = c(method = "kfold", folds = 10).
#'   \item rep_kfold: Random partition in repeated k-fold cross-validation. Usage method = c(method = "rep_kfold", folds = 10, replicates=10). 'folds' refers to the number of folds for data partitioning, it assumes value >=1. 'replicate' refers to the number of replicates, it assumes a value >=1.
#'   \item loocv: Leave-one-out cross-validation (a.k.a. Jackknife) WRITE METHOD HERE . Usage method = c(method = "loocv").
#'   \item boot: Random bootstrap partition. Usage part=c(method='boot', replicates='2', proportion='0.7'). 'replicate' refers to the number of replicates, it assumes a value >=1. 'proportion' refers to the proportion of occurrences used for model fitting, it assumes a value >0 and <=1. In this example proportion='0.7' mean that 70% of data will be used for model training, while 30% for model testing. For this method function will returns .partX columns with "train" and "test" words
#'   }
#'
#' @return
#' A tibble object with information used in 'data' and/or 'gb_data' arguments and a additional columns .part with partition group. In the case of boot partition partition data will be stored with 'train' and 'test' words, for all other methods partition groups will be stored as numeric. In case of use 'data' and 'bg_data' will returned a list with partition data for each database with the names data and bg_data.
#' @export
#'
#' @importFrom dplyr %>% group_by mutate n summarise filter select slice_sample full_join left_join tibble
#'
#' @seealso \code{\link{part_sblock}}, \code{\link{part_senv}}, \code{\link{sample_pseudoabs}}, \code{\link{sample_background}}
#'
#' @examples
#' \dontrun{
#' data("abies_db")
#' abies_db$partition <- NULL
#' abies_db <- tibble(abies_db)
#'
#' # K-fold method
#' abies_db2 <- part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "kfold", folds = 10)
#' )
#' abies_db2
#'
#' # Repeated K-fold method
#' abies_db2 <- part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "rep_kfold", folds = 10, replicates = 10)
#' )
#' abies_db2
#'
#' # Leave-one-out cross-validation (loocv) method
#' abies_db2 <- part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "loocv")
#' )
#' abies_db2
#'
#' # Bootstrap method
#' abies_db2 <- part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "boot", replicates = 50, proportion = 0.7)
#' )
#' abies_db2
#' abies_db2$.part1 %>% table() # Note that for this method .partX columns have train and test words.
#' }
#'
part <- function(data, pr_ab, bg_data = NULL, bg_a = NULL, method = NULL) {
  # TODO add conditional for testing misuse of aguments
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
    filt <- data$.part > filtmi$max
    data$.part[filt] <- rep(1:filtmi$max, length.out = sum(filt))
  }

  # BOOOTSRAP
  if (method["method"] == "boot") {
    reps <- as.numeric(method["replicates"])
    prop <- as.numeric(method["proportion"])
    prop2 <- 1 - prop
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
