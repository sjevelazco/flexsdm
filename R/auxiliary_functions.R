##%######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
##%######################################################%##


pre_tr_te <- function(data, p_names, h){
  train <- list()
  test <- list()

  if(any(c("train", "train-test", "test")
         %in%
         unique(data[, p_names[h]]))){
    np2 <- 1

    filt <- grepl('train', data[,p_names[h]])
    train[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])

    filt <- grepl('test', data[,p_names[h]])
    test[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])

  } else {
    np2 <- max(data[p_names[h]])

    for (i in 1:np2) {

      train[[i]] <- data[data[p_names[h]] == i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])

      test[[i]] <- data[data[p_names[h]] != i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
    }
  }
  return(list(train=train, test=test, np2=np2))
}
