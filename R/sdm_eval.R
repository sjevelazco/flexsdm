#' Calculate different model performance metrics
#'
#' @description This function calculates threshold dependent and independent performance metric. It calculates TPR (True Positive Rate, also called sensitivity), TNR (True Negative Rate, , also called specificity), SORENSEN, JACCARD, FPB, OR (Omission Rate), TSS (True Skill Statistic) KAPPA, AUC (Area Under Curve), BOYCE, IMAE (Inverse Mean Absolute Error, i.e. 1-Mean Absolute Error).
#' @param p numeric. Predicted suitability for presences
#' @param a numeric. Predicted suitability for presences absences
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   }
#' In the case of use more than one threshold type it is necessary concatenate threshold types, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified
#'
#' @param bg numeric. Predicted suitability for background points. It is used for BOYCE metric. It bg is set as null BOYCE metric will be calculated with presences and absences suitability values
#'
#' @return a list with two tibbles
#' \itemize{
#' \item "all_values". A tibble object with the performance metric for 90 threshold values >0 and <1.
#' \item "all_thresholds". A tibble with the threshold values and the performance metric values for all thresholds.
#' \item "selected_thresholds". It is similar to "all_threshold" with difference that it contains the thresholds and the performance metric values for the selected threshold.
#' }
#' when thr argument is NULL function will return a list with "all_values", "all_thresholds" tibbles, when thr argument is used with one or more threshold function returns "selected_thresholds" and "all_thresholds" tibbles
#' By default all thresholds will be calculated
#' @export
#'
#' @importFrom dplyr tibble mutate filter pull bind_cols left_join all_of
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
#' # I the case it is needed use background for calculate all other metric background values can be used in "a" argument
#' sdm_eval(p, backg, thr = "max_sens_spec")
#' }
#'
sdm_eval <- function(p, a, bg = NULL, thr = NULL) {
  if (any(
    !(thr[names(thr) != "sens"]) %in% c(
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

  if (any(thr[grep("type", names(thr))] %in% "sensitivity") &&
    !any(names(thr) %in% "sens")) {
    # stop(
    #   "provide a sensitivity value in the vector used in 'thr' argument, e.g. thr=c(type=c('lpt', 'max_sens_spec', 'sensitivity'), sens='0.8')"
    # )
    thr <- c(thr, sens = 0.9)
  }
  if (!any(thr[grep("type", names(thr))] %in% "sensitivity")) {
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
  if (any(thr[grep("type", names(thr))] %in% "sensitivity")) {
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
    performance <- performance %>% dplyr::mutate(BOYCE = dplyr::mutate(pres = p, contrast = c(p, a)))
  } else {
    performance <- performance %>% dplyr::mutate(BOYCE = dplyr::mutate(pres = p, contrast = c(p, bg)))
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

  thresholds$lpt <- max(performance$threshold[performance$TPR == 1])

  thresholds$sensitivity <- performance$threshold[which(abs(
    performance$TPR - as.numeric(thr["sens"])
  ) ==
    min(abs(performance$TPR - as.numeric(thr["sens"]))))] %>% max()

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
