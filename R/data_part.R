#' Data partitioning for training and testing models
#'
#' @description This function provides different partition methods based in folds (KFOLD, REP_KFOLD, and LOOCV), bootstrap (BOOT), and spatially structured (BANDS and BLOCK).
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records for a given species
#' @param p_a character. Column name of "data" with presences, presence-absence, or pseudo-absence. Presences must be represented by 1 and absences by 0
#' @param bg_data data.frame. Data frames with background points.
#' @param bg_a character. Column name of "bg_data" with absences. Background must be repented by 0
#' @param method character. Vector with partition method to be used:
#' \itemize{
#'   \item KFOLD: Random partition in k-fold cross-validation. Usage part=c(method= 'KFOLD', folds='5'). 'folds' refers to the number of folds for data partitioning, it assumes value >=1. Usage method = c(method = "KFOLD", folds = 10).
#'   \item REP_KFOLD: Random partition in repeated k-fold cross-validation. Usage method = c(method = "REP_KFOLD", folds = 10, replicates=10). 'folds' refers to the number of folds for data partitioning, it assumes value >=1. 'replicate' refers to the number of replicates, it assumes a value >=1.
#'   \item LOOCV: Leave-one-out cross-validation (a.k.a. Jackknife) WRITE METHOD HERE . Usage method = c(method = "LOOCV").
#'   \item BOOT: Random bootstrap partition. Usage part=c(method='BOOT', replicates='2', proportion='0.7'). 'replicate' refers to the number of replicates, it assumes a value >=1. 'proportion' refers to the proportion of occurrences used for model fitting, it assumes a value >0 and <=1. In this example proportion='0.7' mean that 70% of data will be used for model training, while 30% for model testing. For this method function will returns .partX columns with "train" and "test" words
#'   \item BANDS: Geographic partition structured as bands arranged in a latitudinal way (type 1) or longitudinal way (type 2). Usage part=c(method= 'BANDS', type='1'). 'type' refers to the bands disposition
#'   \item BLOCK: Geographic partition structured as a checkerboard (a.k.a. block cross-validation). Usage part=c(method= 'BLOCK').
#'   }
#'
#' @return
#' This function will return the same data.frame used in the arguments 'data' and 'gb_data' with the additional column or columns staring with the names .part with partition group. In the case of BOOT partition data will be stored with 'train' and 'test' words, for all other methods partition groups will be stored as numeric. In case of use 'data' and 'bg_data' will returned a list with partition data for each database with the names data and bg_data.
#' @export
#'
#' @importFrom dplyr group_by mutate n summarise filter select slice_sample full_join left_join tibble
#'
#' @seealso \code{\link{pseudoabs}}, \code{\link{backgroudp}}}
#'
#' @examples
#' \dontrun{
#' require(devtools)
#' data("abies_db")
#' abies_db$partition <- NULL
#' abies_db <- tibble(abies_db)
#'
#' # K-fold method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "KFOLD", folds = 10)
#' )
#' abies_db2
#'
#' # Repeated K-fold method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   bg_data = NULL,
#'   bg_p_a = NULL,
#'   method = c(method = "REP_KFOLD", folds = 10, replicates = 10)
#' )
#' abies_db2
#'
#' # Leave-one-out cross-validation (LOOCV) method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "LOOCV")
#' )
#' abies_db2
#'
#' # Bootstrap method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "BOOT", replicates = 50, proportion = 0.7)
#' )
#' abies_db2$.part1 %>% table() # Note that for this method .partX columns have train and test words.
#' }
#'
data_part <- function(data, p_a, bg_data = NULL, bg_a = NULL, method = NULL) {

  # KFOLD
  if (method["method"] == "KFOLD") {
    data <- data %>%
      dplyr::group_by(!!as.symbol(p_a)) %>%
      dplyr::mutate(.part = sample(rep(1:method["folds"], length.out = dplyr::n()))) %>%
      dplyr::group_by()
  }

  # REP_KFOLD
  if (method["method"] == "REP_KFOLD") {
    for (i in 1:method["replicates"]) {
      cname <- paste0(".part", i)
      data <-
        data %>%
        dplyr::group_by(!!as.symbol(p_a)) %>%
        dplyr::mutate(REP_KFOLD = sample(rep(1:method["folds"],
          length.out = dplyr::n()
        ))) %>%
        dplyr::group_by()
      colnames(data)[colnames(data) == "REP_KFOLD"] <- cname
    }
  }

  # LOOCV
  if (method["method"] == "LOOCV") {
    data <- data %>%
      dplyr::group_by(!!as.symbol(p_a)) %>%
      dplyr::mutate(.part = 1:dplyr::n()) %>%
      dplyr::group_by()
    filt <- data %>%
      dplyr::group_by(!!as.symbol(p_a)) %>%
      dplyr::summarise(max = max(.part))
    filtmi <- filt %>% dplyr::filter(max == min(max))
    filt <- data$.part > filtmi$max
    data$.part[filt] <- rep(1:filtmi$max, length.out = sum(filt))
  }

  # BOOOTSRAP
  if (method["method"] == "BOOT") {
    reps <- as.numeric(method["replicates"])
    prop <- as.numeric(method["proportion"])
    prop2 <- 1 - prop
    data <- data %>% dplyr::mutate(IDBOOT = 1:nrow(data))
    for (i in 1:reps) {
      data2 <- data %>% dplyr::select({{ p_a }}, "IDBOOT")
      data_train <- data2 %>%
        dplyr::group_by(!!as.symbol(p_a)) %>%
        dplyr::slice_sample(prop = prop) %>%
        dplyr::mutate(BOOT1 = "train")
      data_ttest <- data2 %>%
        dplyr::group_by(!!as.symbol(p_a)) %>%
        dplyr::slice_sample(prop = prop2) %>%
        dplyr::mutate(BOOT2 = "test")

      data2 <- dplyr::full_join(data_train, data_ttest, by = c("pr_ab", "IDBOOT")) %>%
        dplyr::mutate(BOOT = paste(BOOT1, BOOT2, sep = "-")) %>%
        dplyr::select(-BOOT1, -BOOT2) %>%
        dplyr::mutate(BOOT = gsub(paste(c("-NA", "NA-"), collapse = "|"), "", BOOT))

      data <- dplyr::left_join(data, data2, by = c("pr_ab", "IDBOOT"))
      colnames(data)[colnames(data) == "BOOT"] <- paste0(".part", i)
    }
    data$IDBOOT <- NULL
  }

  return(dplyr::tibble(data))
}
