#' Title
#'
#' @param pres 
#' @param contrast 
#' @param n_bins 
#' @param n_width 
#' @description This function calculate Boyce index performance metric. Codes were adapted from enmSdm package. Boyce have value between -1 and +1, 
#' with a value tending toward +1 indicating good to perfect predictions, values
#' around 0 indicating predictions no different from those obtained by chance,
#' and values toward -1 indicating counter-predictions. 
#' @return
#' @export
#'
#' @examples
boyce <- function(pres,
                  contrast,
                  n_bins = 101,
                  n_width = 0.1) {
  lowest <- min(c(pres, contrast), na.rm = TRUE)
  highest <-  max(c(pres, contrast), na.rm = TRUE) + .Machine$double.eps
  window_width <- n_width * (highest - lowest)
  
  lows <- seq(lowest, highest - window_width, length.out = n_bins)
  highs <- seq(lowest + window_width + .Machine$double.eps, highest, length.out = n_bins)
  
  ## initiate variables to store predicted/expected (P/E) values
  freq_pres <-  NA
  freq_contrast <- NA
  
  ### tally proportion of test presences/background sites in each class
  for (i in 1:n_bins) {
    # number of presence predictions in a class
    freq_pres[i] <-
      sum(pres >= lows[i] & pres < highs[i], na.rm = TRUE)
    
    # number of background predictions in this class
    freq_contrast[i] <-
      sum(contrast >= lows[i] & contrast < highs[i], na.rm = TRUE)
    
  }
  
  # mean bin prediction
  mean_pred <- rowMeans(cbind(lows, highs))
  
  # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
  if (any(freq_pres > 0 & freq_contrast == 0)) {
    small_value <- 0.5
    freq_contrast[freq_pres > 0 & freq_contrast == 0] <- small_value
  }
  
  # remove classes with 0 presence frequency
  if (any(freq_pres == 0)) {
    zeros <- which(freq_pres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }
  
  # remove classes with 0 background frequency
  if (any(0 %in% freq_contrast)) {
    zeros <- which(freqPres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }
  
  P <- freq_pres / length(pres)
  E <- freq_contrast / length(contrast)
  PE <- P / E
  
  # remove NAs
  rm_nas <- complete.cases(data.frame(mean_pred, PE))
  mean_pred <- mean_pred[rm_nas]
  PE <- PE[rm_nas]
  
  # calculate Boyce index
  result <- stats::cor(x = mean_pred, y = PE, method = 'spearman')
  return(result)
}
