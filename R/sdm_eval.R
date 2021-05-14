#' Calculate different model performance metrics
#'
#' @description This function calculates threshold dependent and independent performance metric. It calculates TPR (True Positive Rate, also called sensitivity), TNR (True Negative Rate, , also called specificity), SORENSEN, JACCARD, FPB, OR (Omission Rate), TSS (True Skill Statistic) KAPPA, AUC (Area Under Curve), BOYCE, IMAE (Inverse Mean Absolute Error, i.e. 1-Mean Absolute Error).
#' @param p numeric. Predicted suitability for presences
#' @param a numeric. Predicted suitability for presences absences
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission. Usage thr=c(type='lpt').
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal (aka threshold that maximizes the TSS).
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest.
#'   Usage thr=c(type='max_sens_spec').
#'   \item max_kappa: The threshold at which Kappa is the highest ("max kappa"). Usage thr=c(type='max_kappa').
#'   \item max_jaccard: The threshold at which Jaccard is the highest. Usage thr=c(type='max_jaccard').
#'   \item max_sorensen: The threshold at which Sorensen is highest. Usage thr=c(type='max_sorensen').
#'   \item max_fpb: The threshold at which FPB is highest. Usage thr=c(type='max_fpb').
#'   \item specific: A threshold value specified by user. Usage thr=c(type='specific', sens='0.6'). 'sens' refers to models will be binarized using this suitability value.
#'   }
#' In the case of use more than one threshold type it is necessary concatenate the names of threshold types, e.g., thr=c(type=c('lpt', 'max_sens_spec', 'max_jaccard')). When specific threshold is used in combination with other it is necessary specify the desired sensitivity value, e.g. thr=c(type=c('lpt', 'max_sens_spec', 'specific'), sens='0.8')
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
#' @export
#'
#' @importFrom dismo evaluate threshold
#' @importFrom dplyr bind_cols left_join
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
#' # Use function without threshold specification
#' e <- sdm_eval(p, a)
#' e$all_values
#' e$threshold
#' e$threshold_table
#'
#' sdm_eval(p, a, thr = c(type = c("max_kappa")))
#' sdm_eval(p, a, thr = c(type = c("lpt", "max_sens_spec", "max_jaccard")))
#' # sdm_eval(p, a, thr = c(type = c("lpt", "max_sens_spec", "specific"))) # wrong way to specific threshold
#' sdm_eval(p, a, thr = c(type = c("lpt", "max_sens_spec", "specific"), sens = "0.8")) # correct way to use specific threshold
#'
#' # Use of bg argument (it will only be used for calculating BOYCE index)
#' sdm_eval(p, a, thr = c(type = c("max_sens_spec")))[[1]]
#' sdm_eval(p, a, thr = c(type = c("max_sens_spec")), bg = backg)[[1]]
#' # I the case it is needed use background for calculate all other metric background values can be used in "a" argument
#' sdm_eval(p, backg, thr = c(type = c("max_sens_spec")))[[1]]
#' }
#'
sdm_eval <- function(p, a, bg = NULL, thr = NULL) {
  if (any(
    !thr[grep("type", names(thr))] %in% c(
      "lpt",
      "max_sens_spec",
      "max_kappa",
      "equal_sens_spec",
      "specific",
      "max_jaccard",
      "max_sorensen",
      "max_fpb"
    )
  )) {
    stop("'thr' Argument is not valid!")
  }

  if (any(thr[grep("type", names(thr))] %in% "specific") &&
    !any(names(thr) %in% "sens")) {
    stop(
      "provide a sensitivity value in the vector used in 'thr' argument, e.g. thr=c(type=c('lpt', 'max_sens_spec', 'specific'), sens='0.8')"
    )
  }

  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop("cannot evaluate a model without absence and presence data that are not NA")
  }

  # Threshold breaks:
  # if(is.null(thr)){

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
  if (any(thr[grep("type", names(thr))] %in% "specific")) {
    tr <- c(tr, as.numeric(thr["sens"]))
    tr <- sort(tr)
  }
  eval_dismo <- dismo::evaluate(p = p, a = a, tr = tr)

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

  # Sorensen Index
  SOR <- 2 * res$tp / (res$fn + (2 * res$tp) + res$fp)
  SorTHR <- tr[which(SOR == max(SOR))][1]

  # Jaccard Index
  JAC <- res$tp / (res$fn + res$tp + res$fp)
  JacTHR <- tr[which(JAC == max(JAC))][1]

  # FPB
  FPB <- 2 * JAC
  FpbTHR <- tr[which(FPB == max(FPB))][1]

  # TSS
  TPR <- res$tp / (res$tp + res$fn)
  TNR <- res$tn / (res$tn + res$fp)
  TSS <- (TPR + TNR) - 1
  TSSTHR <- tr[which(TSS == max(TSS))][1]

  # Threshold
  thresholds <- list()
  thresholds$max_sorensen <- SorTHR
  thresholds$max_jaccard <- JacTHR
  thresholds$max_fpb <- JacTHR

  ThrDis <- c("kappa", "spec_sens", "no_omission", "equal_sens_spec", "sensitivity")
  ThrDis <- (sapply(ThrDis, function(x) {
    dismo::threshold(eval_dismo)[x]
  }))
  names(ThrDis) <- nom <- c("max_kappa", "max_sens_spec", "lpt", "equal_sens_spec", "specific")
  if (any(thr[grep("type", names(thr))] %in% "specific")) {
    ThrDis$specific <- as.numeric(thr["sens"])
  }
  thresholds <- c(ThrDis, thresholds)
  thresholds <- dplyr::bind_cols(thresholds)

  # Final Result Object
  performance <- list()
  performance$threshold <- tr
  performance$n_presences <- np
  performance$n_absences <- na
  performance$TPR <- res$tp / (res$tp + res$fn)
  performance$TNR <- res$tn / (res$tn + res$fp)
  performance$SORENSEN <- SOR
  performance$JACCARD <- JAC
  performance$FPB <- FPB
  performance$OR <- (1 - performance$TPR)
  performance$TSS <- TSS
  performance$KAPPA <- eval_dismo@kappa
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1) / 2)
  performance$AUC <- R / (as.numeric(na) * as.numeric(np))

  if (is.null(bg)) {
    performance$BOYCE <- boyce(pres = p, contrast = c(p, a))
  } else {
    performance$BOYCE <- boyce(pres = p, contrast = c(p, bg))
  }
  real <- c(rep(1, length(p)), rep(0, length(a)))
  pred <- c(p, a)
  performance$IMAE <- 1 - (sum(abs(real - pred)) / length(pred))

  performance <- dplyr::bind_cols(performance)

  thr_table <- dplyr::tibble(threshold = names(thresholds), values = unlist(thresholds))
  thr_table <- dplyr::left_join(thr_table, performance, by = c("values" = "threshold"))

  # Return final result

  result <-
    list(
      all_values = performance,
      all_thresholds = thr_table
    )

  if (is.null(thr)) {
    return(result)
  } else {
    result <- result$all_threshold
    result <-
      list(
        selected_thresholds = result[result$threshold %in% thr, ],
        all_thresholds = result
      )
    return(result)
  }
}
