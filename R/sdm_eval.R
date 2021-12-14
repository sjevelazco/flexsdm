#' Calculate different model performance metrics
#'
#' @description This function calculates threshold dependent and independent model performance
#' metrics.
#'
#' @param p numeric. Predicted suitability for presences
#' @param a numeric. Predicted suitability for absences
#' @param thr character. Threshold criterion used to get binary suitability values (i.e. 0,1).
#' Used for threshold-dependent performance metrics.
#' It is possible to use more than one threshold type.
#' A vector must be provided for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the Sensitivity and Specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the Sensitivity and Specificity
#'   is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is the highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is the highest.
#'   \item sensitivity: Threshold based on a specified Sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers
#'   to Sensitivity value. If a sensitivity value is not specified, the
#'    default value is 0.9
#'   }
#' If more than one threshold type is used, concatenate threshold types,
#' e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec',
#' 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity').
#' Function will use all thresholds if no threshold type is specified
#' @param bg numeric. Predicted suitability for background points, used for BOYCE metric.
#' It bg is set as NULL, BOYCE metric will be calculated with presences and absences suitabilities
#' values
#'
#' @details This function is used for evaluating different models approaches base on the combination
#' of presence-absences or presence-pseudo-absences and background point data and
#' suitability predicted by any model or flexsdm modeling function families (fit_, esm_, and tune_.)
#'
#' It calculates the next performance metric:
#'
#'   | Performance metric | Threshold dependent   | Values ranges  |
#'   | :------------- |:-------------:| -----:|
#'   | TPR (True Positive Rate, also called Sensitivity) | yes | 0 - 1 |
#'   | TNR (True Negative Rate, also called Specificity) | yes | 0 - 1 |
#'   | SORENSEN                                          | yes | 0 - 1 |
#'   | JACCARD                                           | yes | 0 - 1 |
#'   | FPB (F-measure on presence-background)            | yes | 0 - 2 |
#'   | OR (Omission Rate)                                | yes | 0 - 1 |
#'   | TSS (True Skill Statistic)                        | yes | -1 - 1 |
#'   | KAPPA                                             | yes | 0 - 1 |
#'   | AUC (Area Under Curve)                            | no | 0 - 1 |
#'   | BOYCE  (continuous Boyce index)*                  | no | -1 - 1 |
#'   | IMAE (Inverse Mean Absolute Error)**              | no | 0 - 1 |
#'
#' \* BOYCE is calculated based on presences and background points, in case that background points
#' is not provided it is calculated using presences and absences. The codes for calculating
#' this metric is and adaptation of enmSdm package (https://github.com/adamlilith/enmSdm)
#'
#' \** IMAE is calculated as 1-(Mean Absolute Error) in order to be consistent with the other
#' metrics where the higher the value of a given performance metric, the greater the model's
#' accuracy
#'
#' @md
#'
#' @return a tibble with next columns
#' \itemize{
#' \item threshold: threshold names
#' \item thr_value: threshold values
#' \item n_presences: number of presences
#' \item n_absences: number of absences
#' \item from TPR to IMAE: performance metrics
#' }
#'
#'
#' @export
#'
#' @importFrom dplyr %>% tibble mutate filter pull bind_cols left_join all_of
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#'
#' set.seed(0)
#' p <- rnorm(50, mean = 0.7, sd = 0.3) %>% abs()
#' p[p > 1] <- 1
#' p[p < 0] <- 0
#'
#' set.seed(0)
#' a <- rnorm(50, mean = 0.3, sd = 0.2) %>% abs()
#' a[a > 1] <- 1
#' a[a < 0] <- 0
#'
#' set.seed(0)
#' backg <- rnorm(1000, mean = 0.4, sd = 0.4) %>% abs()
#' backg[backg > 1] <- 1
#' backg[backg < 0] <- 0
#'
#' # Function use without threshold specification
#' e <- sdm_eval(p, a)
#' e
#'
#' # Function use with threshold specification
#' sdm_eval(p, a, thr = "max_sorensen")
#' sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "max_jaccard"))
#' sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "sensitivity"))
#' sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "sensitivity", sens = "0.95"))
#'
#' # Use of bg argument (it will only be used for calculating BOYCE index)
#' sdm_eval(p, a, thr = "max_sens_spec")
#' sdm_eval(p, a, thr = c("max_sens_spec"), bg = backg)
#'
#' # If background will be used to calculate all other metrics
#' # background values can be used in "a" argument
#' sdm_eval(p, backg, thr = "max_sens_spec")
#' }
#'
sdm_eval <- function(p, a, bg = NULL, thr = NULL) {
  TPR <- TNR <- JACCARD <- SORENSEN <- threshold <- FPB <- TSS <- NULL
  if (any(
    !(thr[is.na(suppressWarnings(as.numeric(thr)))]) %in% c(
      "lpt",
      "max_sens_spec",
      "equal_sens_spec",
      "sensitivity",
      "max_jaccard",
      "max_sorensen",
      "max_fpb"
    )
  )) {
    stop("'thr' Argument is not valid!")
  }

  if (is.null(thr)) {
    thr <- c(
      "lpt",
      "max_sens_spec",
      "max_kappa",
      "equal_sens_spec",
      "sensitivity",
      "max_jaccard",
      "max_sorensen",
      "max_fpb"
    )
  }

  if (any(thr %in% "sensitivity") &&
    !any(names(thr) %in% "sens")) {
    thr <- c(thr, sens = 0.9)
  }

  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop("cannot evaluate a model without absence and presence data that are not NA")
  }

  # Threshold breaks:
  if (length(p) > 1000) {
    tr <- as.vector(stats::quantile(p, 0:1000 / 1000))
  } else {
    tr <- p
  }
  if (length(a) > 1000) {
    tr <- c(tr, as.vector(stats::quantile(a, 0:1000 / 1000)))
  } else {
    tr <- c(tr, a)
  }
  tr <- sort(unique(round(tr, 8)))
  if (any(thr %in% "sensitivity")) {
    tr <- c(tr, as.numeric(thr["sens"]))
    tr <- sort(tr)
  }

  res <- matrix(ncol = 4, nrow = length(tr))
  colnames(res) <- c("tp", "fp", "fn", "tn")
  # Confusion Matrix
  for (i in 1:length(tr)) {
    res[i, 1] <- length(p[p >= tr[i]]) # a  true positives
    res[i, 2] <- length(a[a >= tr[i]]) # b  false positives
    res[i, 3] <- length(p[p < tr[i]]) # c  false negatives
    res[i, 4] <- length(a[a < tr[i]]) # d  true negatives
  }
  res <- data.frame(res)

  # Performance metrics
  performance <- dplyr::tibble(
    threshold = tr,
    n_presences = np,
    n_absences = na,
    TPR = res$tp / (res$tp + res$fn),
    TNR = res$tn / (res$tn + res$fp),
    SORENSEN = 2 * res$tp / (res$fn + (2 * res$tp) + res$fp),
    JACCARD = res$tp / (res$fn + res$tp + res$fp),
    FPB = 2 * JACCARD,
    OR = (1 - TPR),
    TSS = (TPR + TNR) - 1
  )

  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1) / 2)
  performance <- performance %>% dplyr::mutate(AUC = R / (as.numeric(na) * as.numeric(np)))

  if (is.null(bg)) {
    performance <- performance %>% dplyr::mutate(BOYCE = boyce(pres = p, contrast = c(p, a)))
  } else {
    performance <- performance %>% dplyr::mutate(BOYCE = boyce(pres = p, contrast = c(p, bg)))
  }

  real <- c(rep(1, length(p)), rep(0, length(a)))
  pred <- c(p, a)
  performance <- performance %>% dplyr::mutate(IMAE = 1 - (sum(abs(real - pred)) / length(pred)))

  # Thresholds
  thresholds <- list()

  thresholds$max_sorensen <-
    max(performance %>% dplyr::filter(SORENSEN == max(SORENSEN)) %>% dplyr::pull(threshold))

  thresholds$max_jaccard <-
    max(performance %>% dplyr::filter(SORENSEN == max(SORENSEN)) %>% dplyr::pull(threshold))

  thresholds$max_fpb <-
    max(performance %>% dplyr::filter(FPB == max(FPB)) %>% dplyr::pull(threshold))

  thresholds$max_sens_spec <-
    max(performance %>% dplyr::filter(TSS == max(TSS)) %>% dplyr::pull(threshold))

  thresholds$equal_sens_spec <-
    performance$threshold[which(abs(
      performance$TPR - performance$TNR
    ) ==
      min(abs(performance$TPR - performance$TNR)))] %>% max()

  suppressWarnings(thresholds$lpt <- max(performance$threshold[performance$TPR == 1]))

  if (any(thr == "sensitivity")) {
    thresholds$sensitivity <- performance$threshold[which(abs(
      performance$TPR - as.numeric(thr["sens"])
    ) ==
      min(abs(performance$TPR - as.numeric(thr["sens"]))))] %>% max()
  }


  thresholds <- dplyr::bind_cols(thresholds)

  thr_table <- dplyr::tibble(threshold = names(thresholds), thr_value = unlist(thresholds))
  thr_table <- dplyr::left_join(thr_table, performance, by = c("thr_value" = "threshold"))

  # Return final result
  result <- thr_table
  if (!is.null(thr)) {
    result <- result %>%
      dplyr::filter(threshold %in% dplyr::all_of(thr))
  }

  return(result)
}
