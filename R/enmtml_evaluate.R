#' Calculate different model peformance metrics
#'
#' @param p 
#' @param a 
#' @param thr 
#' @param bg 
#'
#' @return
#' @export
#' 
#' @importFrom dismo evaluate threshold
#' @importFrom dplyr bind_cols left_join
#' 
#' @examples
enmtml_evaluate <- function(p, a, bg=NULL, thr=NULL){
  #Parameters:
  #p:presence points suitability
  #a:absence points suitability
  #tr:numeric vector with threshold values
  
  require(dismo)
  require(dplyr)
  
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
  
  #Fpb
  Fpb <- 2 * JAC
  # if(is.null(thr)){
    FpbTHR <- tr[which(Fpb == max(Fpb))][1]
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
  performance$Fpb <- Fpb
  performance$OR  <- (1-performance$TPR)
  performance$TSS <- (performance$TPR + performance$TNR) - 1
  performance$KAPPA <- eval_dismo@kappa
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
  performance$AUC <- R/(as.numeric(na) * as.numeric(np))
  
  if(is.null(bg)){
    performance$BOYCE <- boyce(pres = p, contrast = c(p, a))
  } else {
    performance$BOYCE <- boyce(pres = p, contrast = c(p, bg))
  }
  
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

