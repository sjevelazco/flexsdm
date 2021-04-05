#' Calculate different model performance metrics
#'
#' @description This function calculates threshold dependent and independent performance metric. It calculates TPR (True Positive Rate, also called sensitivity), TNR (True Negative Rate, , also called specificity), SORENSEN, JACCARD, FPB, OR (Omission Rate), TSS (True Skill Statistic) KAPPA, AUC (Area Under Curve), BOYCE, MAE (Mean Absolute Error, values closer to zero denote better performance).
#' @param p numeric. Predicted suitability for presences
#' @param a numeric. Predicted suitability for presences absences
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item LPT: The highest threshold at which there is no omission. Usage thr=c(type='LPT').
#'   \item EQUAL_SENS_SPEC: Threshold at which the sum of the sensitivity and specificity is the highest.
#'   \item MAX_TSS: Threshold at which the sensitivity and specificity are equal.
#'   Usage thr=c(type='MAX_TSS').
#'   \item MAX_KAPPA: The threshold at which kappa is the highest ("max kappa"). Usage thr=c(type='MAX_KAPPA').
#'   \item MAX_JACCARD: The threshold at which Jaccard is the highest. Usage thr=c(type='MAX_JACCARD').
#'   \item MAX_SORENSEN: The threshold at which Sorensen is highest. Usage thr=c(type='MAX_SORENSEN').
#'   \item MAX_FPB: The threshold at which FPB is highest. Usage thr=c(type='MAX_FPB').
#'   \item SENSITIVITY: A threshold value specified by user. Usage thr=c(type='SENSITIVITY', sens='0.6'). 'sens' refers to models will be binarized using this suitability value.
#'   }
#' In the case of use more than one threshold type it is necessary concatenate the names of threshold types, e.g., thr=c(type=c('LPT', 'MAX_TSS', 'MAX_JACCARD')). When SENSITIVITY threshold is used in combination with other it is necessary specify the desired sensitivity value, e.g. thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')
#'
#' @param bg numeric. Predicted suitability for background points. It is used for BOYCE metric. It bg is set as null BOYCE metric will be calculated with presences and absences suitability values
#'
#' @return a list with the next tibble
#' \itemize{
#' \item "performance". A tibble object with the performance metric for 90 threshold values >0 and <1. Performance.
#' \item "threshold".  A tibble with values for all the threshold available in the package.
#' \item "threshold_table". A tibble with the threshold values and the performance metric values for each threshold.
#' \item " selected_threshold". It is similar to "threshold_table" with difference that it contains the thresholds and the performance metric values for the selected threshold. A tibble with threshold values and the performance metric values for each threshold.
#' }
#' when thr argument is NULL function will return a list with "performance", "threshold", and "threshold_table" tibbles, when thr argument is used with one or more threshold function returns "selected_threshold" and "threshold_table" tibbles
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
#' p <- rnorm(50, mean=0.7, sd=0.3) %>% abs()
#' p[p>1] <- 1
#' p[p<0] <- 0
#'
#' set.seed(0)
#' a <- rnorm(50, mean=0.3, sd=0.2) %>% abs()
#' a[a>1] <- 1
#' a[a<0] <- 0
#'
#' set.seed(0)
#' backg <- rnorm(1000, mean=0.4, sd=0.4) %>% abs()
#' backg[backg>1] <- 1
#' backg[backg<0] <- 0
#'
#' # Use function without threshold specification
#' e <- enm_eval(p, a)
#' e$performance
#' e$threshold
#' e$threshold_table
#'
#' enm_eval(p, a, thr=c(type=c('MAX_KAPPA')))
#' enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'MAX_JACCARD')))
#' enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'))) # wrong way to SENSITIVITY threshold
#' enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')) # correct way to use SENSITIVITY threshold
#'
#' # Use of bg argument (it will only be used for calculating BOYCE index)
#' enm_eval(p, a, thr=c(type=c('MAX_TSS')))[[1]]
#' enm_eval(p, a, thr=c(type=c('MAX_TSS')), bg=backg)[[1]]
#' # I the case it is needed use background for calculate all other metric background values can be used in "a" argument
#' enm_eval(p, backg, thr=c(type=c('MAX_TSS')))[[1]]
#' }
#'
enm_eval <- function(p, a, bg=NULL, thr=NULL){

  if(any(
    !thr[grep('type', names(thr))] %in% c(
      "LPT",
      "MAX_TSS",
      "MAX_KAPPA",
      "EQUAL_SENS_SPEC",
      "SENSITIVITY",
      "MAX_JACCARD",
      "MAX_SORENSEN",
      "MAX_FPB"
    )
  )) {
    stop("'thr' Argument is not valid!")
  }

  if (any(thr[grep("type", names(thr))] %in% "SENSITIVITY") &&
      !any(names(thr) %in% "sens")) {
    stop(
      "provide a sensitivity value in the vector used in 'thr' argument, e.g. thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')")
  }

  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }

  #Threshold breaks:
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
    if(any(thr[grep("type", names(thr))] %in% "SENSITIVITY")){
      tr <- c(tr, as.numeric(thr['sens']))
      tr <- sort(tr)
    }
    # tr <- c(tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
    # tr <- abs(sort(tr))
    # tr <- tr[tr<=1]
    eval_dismo <- dismo::evaluate(p=p, a=a, tr = tr)
  # }else{
  #   tr <- thr
  #
  #   eval_dismo <- dismo::evaluate(p=p, a=a, tr = tr)
  # }

  res <- matrix(ncol = 4, nrow = length(tr))
  colnames(res) <- c('tp', 'fp', 'fn', 'tn')
  #Confusion Matrix
  for (i in 1:length(tr)) {
    res[i, 1] <- length(p[p >= tr[i]])  # a  true positives
    res[i, 2] <- length(a[a >= tr[i]])  # b  false positives
    res[i, 3] <- length(p[p < tr[i]])    # c  false negatives
    res[i, 4] <- length(a[a < tr[i]])    # d  true negatives
  }
  res <- data.frame(res)

  #Sorensen Index
  SOR <- 2 * res$tp / (res$fn + (2 * res$tp) + res$fp)
  # if(is.null(thr)){
    SorTHR <- tr[which(SOR == max(SOR))][1]
  # }else{
  #   SorTHR <- tr
  # }

  #Jaccard Index
  JAC <- res$tp / (res$fn + res$tp + res$fp)
  # if(is.null(thr)){
    JacTHR <- tr[which(JAC == max(JAC))][1]
  # }else{
  #   JacTHR <- tr
  # }

  #FPB
  FPB <- 2 * JAC
  # if(is.null(thr)){
    FpbTHR <- tr[which(FPB == max(FPB))][1]
  # }else{
  #   FpbTHR <- tr
  # }

  #TSS
  TPR <- res$tp / (res$tp + res$fn)
  TNR <- res$tn / (res$tn + res$fp)
  TSS <- (TPR + TNR) - 1
  # if(is.null(thr)){
  TSSTHR <- tr[which(TSS == max(TSS))][1]
  # }else{
  #   TSSTHR <- tr
  # }

  # # SEDI
  # Hi <- res$tp / (res$tp + res$fn) #true positive rate
  # Fi <- res$fp / (res$tn + res$fp) #false positive rate
  # Hi <- ifelse(Hi == 0, .Machine$double.eps, Hi)
  # Fi <- ifelse(Fi == 0, .Machine$double.eps, Fi)
  # Hi_minus1 <- 1 - Hi
  # Fi_minus1 <- 1 - Fi
  # Hi_minus1 <- ifelse(Hi_minus1 == 0, .Machine$double.eps, Hi_minus1)
  # Fi_minus1 <- ifelse(Fi_minus1 == 0, .Machine$double.eps, Fi_minus1)
  #
  # SEDI <-
  #   (log(Fi)-log(Hi)-log(Fi_minus1) + log(Hi_minus1)) /
  #   (log(Fi)+log(Hi)+log(Fi_minus1) + log(Hi_minus1))
  #
  # # ORSS
  # ORSS <- (res$tp*res$tn-res$fp*res$fn)/(res$tp*res$tn+res$fp*res$fn)


  thresholds <- list()
  thresholds$MAX_SORENSEN <- SorTHR
  thresholds$MAX_JACCARD <- JacTHR
  thresholds$MAX_FPB <- JacTHR

  ThrDis <- c("kappa", "spec_sens", "no_omission", "equal_sens_spec", "sensitivity")
  ThrDis <- (sapply(ThrDis, function(x)
      dismo::threshold(eval_dismo)[x]))
  names(ThrDis) <- nom <- c("MAX_KAPPA", "MAX_TSS", "LPT", "EQUAL_SENS_SPEC", "SENSITIVITY")
  if(any(thr[grep("type", names(thr))] %in% "SENSITIVITY")){
    ThrDis$SENSITIVITY <- as.numeric(thr['sens'])
  }
  thresholds <- c(ThrDis, thresholds)
  thresholds <- dplyr::bind_cols(thresholds)

  #Final Result Object
  performance <- list()
  performance$threshold <- tr
  performance$n_presences <- np
  performance$n_absences <- na
  performance$TPR <- res$tp / (res$tp + res$fn)
  performance$TNR <- res$tn / (res$tn + res$fp)
  performance$SORENSEN <- SOR
  performance$JACCARD <- JAC
  performance$FPB <- FPB
  performance$OR  <- (1-performance$TPR)
  performance$TSS <- TSS
  performance$KAPPA <- eval_dismo@kappa
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
  performance$AUC <- R/(as.numeric(na) * as.numeric(np))

  if(is.null(bg)){
    performance$BOYCE <- boyce(pres = p, contrast = c(p, a))
  } else {
    performance$BOYCE <- boyce(pres = p, contrast = c(p, bg))
  }
  real <- c(rep(1, length(p)), rep(0, length(a)))
  pred <- c(p, a)
  performance$MAE <- sum(abs(real-pred))/length(pred)

  performance <- dplyr::bind_cols(performance)

  thr_table <- dplyr::tibble(threshold = names(thresholds), values=unlist(thresholds))
  thr_table <- dplyr::left_join(thr_table, performance, by=c('values'='threshold'))

  #Return final result

  result <-
    list(
      performance = performance,
      threshold = thresholds,
      threshold_table = thr_table
    )

  if (is.null(thr)) {
    return(result)
  } else {
    result <- result$threshold_table
    result <-
      list(selected_threshold = result[result$threshold %in% thr,],
           threshold_table = result)
    return(result)
  }
}
