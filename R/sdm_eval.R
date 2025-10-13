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
#' If not provided (NULL), the Boyce index will be calculated using absence data instead.
#' **Note:** Using absence data for the Boyce index is not standard practice and may result in inflated performance values. It is highly recommended to provide background points for a correct calculation.
#' The Boyce index is calculated using the `boyce` function, which is an adaptation of the method implemented in the `ecospat` package.
#'
#' @md
#' @details This function is used for evaluating different models approaches base on the combination
#' of presence-absences or presence-pseudo-absences and background point data and
#' suitability predicted by any model or flexsdm modeling function families (fit_, esm_, and tune_.)
#'
#' ### Performance Metrics Formulas
#'
#' It calculates the next performance metric:
#'
#'   | Performance metric | Threshold dependent   | Values ranges  |
#'   | :------------- |:-------------:| -----:|
#'   | TPR (True Positive Rate, also called Sensitivity) | yes | 0 - 1 |
#'   | TNR (True Negative Rate, also called Specificity) | yes | 0 - 1 |
#'   | W_TPR_TNR (Weighted TPR-TNR; Li et al. 2020)      | yes | 0 - 1 |
#'   | SORENSEN                                          | yes | 0 - 1 |
#'   | JACCARD                                           | yes | 0 - 1 |
#'   | FPB (F-measure on presence-background)            | yes | 0 - 2 |
#'   | OR (Omission Rate)                                | yes | 0 - 1 |
#'   | TSS (True Skill Statistic)                        | yes | -1 - 1 |
#'   | KAPPA                                             | yes | 0 - 1 |
#'   | MCC (Matthews Correlation Coefficient; Matthews 1975)            | yes | -1 - 1 (1 is best)         |
#'   | AUC (Area Under Curve)                            | no | 0 - 1 |
#'   | BOYCE  (continuous Boyce index)*                  | no | -1 - 1 |
#'   | IMAE (Inverse Mean Absolute Error)**              | no | 0 - 1 |
#'   | CRPS (Continuous Ranked Probability Score based on Brier Score, Brier 1950)**        | no  | 0 - 1 |
#'
#' \* The continuous Boyce index is calculated based on presences and background points. If background points are not provided, it will be calculated using presences and absences, which is not standard and may lead to misleading results. The code for calculating this metric is an adaptation of the `ecospat` package.
#'
#' \** IMAE and CRPS are calculated as 1-(Mean Absolute Error) and 1-(CRPS), respectively, in order to be consistent with the other
#' metrics where the higher the value of a given performance metric, the greater the model's.
#'
#'
#' To define the formulas, the following components of the confusion matrix are used:
#' - `tp`: True Positives (presences correctly predicted as presences)
#' - `tn`: True Negatives (absences correctly predicted as absences)
#' - `fp`: False Positives (absences incorrectly predicted as presences)
#' - `fn`: False Negatives (presences incorrectly predicted as absences)
#' - `np`: Number of presences (`length(p)`)
#' - `na`: Number of absences (`length(a)`)
#'
#' The formulas are:
#' - **TPR (Sensitivity)**: \deqn{TPR = \frac{tp}{tp + fn}}{TPR = tp / (tp + fn)}
#' - **TNR (Specificity)**: \deqn{TNR = \frac{tn}{tn + fp}}{TNR = tn / (tn + fp)}
#' - **W_TPR_TNR** (Li et al. 2020): \deqn{W\_TPR\_TNR = w \cdot TPR + (1-w) \cdot TNR}{W_TPR_TNR = w * TPR + (1-w) * TNR}, where \deqn{w = \frac{na}{na + np}}{w = na / (na + np)}
#' - **SORENSEN**: \deqn{Sorensen = \frac{2 \cdot tp}{fn + 2 \cdot tp + fp}}{Sorensen = 2*tp / (fn + 2*tp + fp)}
#' - **JACCARD**: \deqn{Jaccard = \frac{tp}{fn + tp + fp}}{Jaccard = tp / (fn + tp + fp)}
#' - **FPB**: \deqn{FPB = 2 \cdot Jaccard}{FPB = 2 * Jaccard}
#' - **OR**: \deqn{OR = 1 - TPR}{OR = 1 - TPR}
#' - **TSS**: \deqn{TSS = TPR + TNR - 1}{TSS = TPR + TNR - 1}
#' - **KAPPA**: \deqn{KAPPA = \frac{Pr(a) - Pr(e)}{1 - Pr(e)}}{KAPPA = (Pr(a) - Pr(e)) / (1 - Pr(e))}, where \deqn{Pr(a) = \frac{tp+tn}{tp+tn+fp+fn}}{Pr(a) = (tp+tn)/(tp+tn+fp+fn)} and \deqn{Pr(e) = \frac{(tp+fp)(tp+fn) + (fn+tn)(fp+tn)}{(tp+tn+fp+fn)^2}}{Pr(e) = ((tp+fp)*(tp+fn) + (fn+tn)*(fp+tn))/(tp+tn+fp+fn)^2}
#' - **MCC** (Matthews 1975): \deqn{MCC = \frac{(tp \cdot tn) - (fp \cdot fn)}{\sqrt{(tp+fp)(tp+fn)(tn+fp)(tn+fn)}}}{MCC = (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))}
#' - **AUC**: Calculated as the Wilcoxon-Mann-Whitney U statistic, which is equivalent to the area under the ROC curve.
#' - **BOYCE**: The continuous Boyce index, which measures how model predictions differ from a random distribution of observed presences across the prediction gradient.
#' - **CRPS** (Brier 1950): For binary outcomes, this is calculated as \deqn{1 - \frac{\sum(predicted - observed)^2}{N}}{1 - (sum(predicted - observed)^2/N)}, which is 1 minus the Brier Score.
#' - **IMAE**: \deqn{IMAE = 1 - \frac{\sum|predicted - observed|}{N}}{IMAE = 1 - (sum|predicted - observed|/N)}, where N is the total number of records.
#'
#' @references
#' \itemize{
#'   \item Brier GW. (1950) Verification of forecasts expressed in terms of probability. Monthly Weather Review 78(1): 1–3. https://doi.org/10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2
#'   \item Li, J., Liu, H., & Li, L. (2020). A novel performance metric for imbalanced learning and its application in credit default prediction. Expert Systems with Applications, 152, 113382. https://doi.org/10.1016/j.eswa.2020.113382
#'   \item Matthews BW. (1975) Comparison of the predicted and observed secondary structure of T4 phage lysozyme. Biochim Biophys Acta (BBA) Protein Struct. 405(2):442–51. https://doi.org/10.1016/0005-2795(75)90109-9
#' }
#'
#' @return a tibble with next columns
#' \itemize{
#' \item threshold: threshold names
#' \item thr_value: threshold values
#' \item n_presences: number of presences
#' \item n_absences: number of absences
#' \item from TPR to CRPS: performance metrics
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
sdm_eval <- function(p, a, bg = NULL, thr = NULL) {
  TPR <- TNR <- JACCARD <- SORENSEN <- threshold <- FPB <- TSS <- NULL

  # Boyce Index based
  # This function calculate Boyce index performance metric.
  # Codes were adapted from ecospat package.
  # Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A. (2006).
  # Evaluating the ability of habitat suitability models to predict species presences.
  # Ecological Modelling, 199(2), 142-152.
  #
  boyce_ <- function(fit, obs, n_bins = 101) {
    # Range of suitability values
    fit.all <- c(fit, obs)
    if (length(unique(fit.all)) == 1) {
      return(NA)
    }

    # Determine window width
    h_w <- (max(fit.all) - min(fit.all)) / 10

    # Create bins
    if (length(unique(fit.all)) < n_bins) {
      bins <- sort(unique(fit.all))
    } else {
      bins <- seq(min(fit.all), max(fit.all), length.out = n_bins)
    }

    # Calculate frequencies
    obs_freq <- sapply(bins, function(i) {
      sum(obs >= i - h_w & obs < i + h_w)
    })

    fit_freq <- sapply(bins, function(i) {
      sum(fit >= i - h_w & fit < i + h_w)
    })

    # Calculate Predicted/Expected ratio
    P <- obs_freq / sum(obs_freq)
    E <- fit_freq / sum(fit_freq)

    # Handle cases where E is zero
    PE_ratio <- P / E

    # Remove bins with no presences
    to_keep <- which(P > 0 & E > 0)

    if (length(to_keep) < 2) {
      return(NA)
    }

    # Calculate Spearman correlation
    cor(bins[to_keep], PE_ratio[to_keep], method = "spearman")
  }

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
    stop("cannot evaluate a model without absence and presence data")
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
  tr <- sort(unique(tr))

  res <- matrix(ncol = 4, nrow = length(tr))
  colnames(res) <- c("tp", "fp", "fn", "tn")
  # Confusion Matrix
  for (i in 1:length(tr)) {
    res[i, 1] <- sum(p >= tr[i]) # a  true positives
    res[i, 2] <- sum(a >= tr[i]) # b  false positives
    res[i, 3] <- sum(p < tr[i]) # c  false negatives
    res[i, 4] <- sum(a < tr[i]) # d  true negatives
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
    KAPPA = (res$tp + res$tn) / (res$tp + res$tn + res$fp + res$fn) -
      (res$fp + res$fn) / (res$tp + res$tn + res$fp + res$fn)^2,
    FPB = 2 * JACCARD,
    OR = (1 - TPR),
    TSS = (TPR + TNR) - 1
  )

  # Add Weighted TPR-TNR
  # w = n / (n + p) -> na / (na + np)
  w <- na / (na + np)
  performance <- performance %>%
    dplyr::mutate(W_TPR_TNR = w * TPR + (1 - w) * TNR)
  # dplyr::mutate(W_TPR_TNR = TPR * (na/(na+np))  + TNR * (np/(na+np)))

  # Add MCC (Matthews Correlation Coefficient)
  mcc_num <- (res$tp * res$tn) - (res$fp * res$fn)
  mcc_den <- sqrt(as.numeric(res$tp + res$fp) * as.numeric(res$tp + res$fn) * as.numeric(res$tn + res$fp) * as.numeric(res$tn + res$fn))
  mcc <- mcc_num / mcc_den
  mcc[is.na(mcc)] <- 0 # Handle cases where the denominator is 0

  performance <- performance %>% dplyr::mutate(MCC = mcc)


  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1) / 2)
  performance <- performance %>% dplyr::mutate(AUC = R / (as.numeric(na) * as.numeric(np)))

  if (is.null(bg)) {
    performance <- performance %>% dplyr::mutate(BOYCE = boyce_(fit = a, obs = p, n_bins = 101))
  } else {
    performance <- performance %>% dplyr::mutate(BOYCE = boyce_(fit = bg, obs = p, n_bins = 101))
  }

  real <- c(rep(1, length(p)), rep(0, length(a)))
  pred <- c(p, a)
  performance <- performance %>% dplyr::mutate(IMAE = 1 - (sum(abs(real - pred)) / length(pred)))

  # Add 1-CRPS (Brier Score for binary outcomes)
  crps_val <- 1 - mean((pred - real)^2)
  performance <- performance %>% dplyr::mutate(CRPS = crps_val)


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

  # finding the maximum threshold
  # that results in a sensitivity >= target
  if (any(thr == "sensitivity")) {
    thresholds$sensitivity <- max(performance$threshold[performance$TPR >=
      as.numeric(thr["sens"])])
  }

  thresholds <- dplyr::bind_cols(thresholds)

  thr_table <- dplyr::tibble(threshold = names(thresholds), thr_value = unlist(thresholds))
  thr_table <- dplyr::left_join(thr_table, performance, by = c("thr_value" = "threshold"))

  # Return final result
  result <- thr_table
  if (!is.null(thr)) {
    result <- result %>%
      dplyr::filter(threshold %in% thr)
  }

  result <- result[c(
    "threshold", "thr_value", "n_presences", "n_absences", "TPR", "TNR", "W_TPR_TNR",
    "SORENSEN", "JACCARD", "FPB", "OR", "TSS", "KAPPA", "MCC", "AUC", "BOYCE", "CRPS", "IMAE"
  )]
  return(result)
}
