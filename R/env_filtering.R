#' Function to perform environmental filtering of species occurrence records
#'
#' @param coord data.frame. A data.frame with longitude and latitude in the first and second columns.
#' @param variables data.frame. A data.frame with environmental conditions. It is possible use two or three variables.  
#' @param nbins numeric. A number of classes used to split each environmental condition.
#' @param plot logical. Plot filtering procedure.
env_filtering <- function(coord, variables, nbins, plot=TRUE){
  library(maps)
  
  filters <- as.list(variables)
  
  res  <- lapply(filters, function(x) {
    (max(x) - min(x)) / nbins
  })
  
  n <- length(filters)
  pot_p <- list ()
  for (i in 1:n) {
    k <- filters [[i]] [!is.na(filters[[i]])]
    ext1 <- range (k)
    ext1[1] <- ext1[1] - 1
    x <- seq(ext1[1], ext1[2], by = res[[i]])
    pot_p[[i]] <- x
  }
  pot_p <- expand.grid(pot_p)
  
  ends <- NULL
  for (i in 1:n) {
    f <- pot_p[, i] + res[[i]]
    ends <- cbind (ends, f)
  }
  
  pot_pp <- data.frame(pot_p, ends)
  pot_pp <- data.frame(pot_pp, groupID = c(1:nrow(pot_pp)))
  rows <- length(filters[[1]])
  # filter <- data.frame(matrix(unlist(filters), nrow = rows))
  real_p <- data.frame(coord, variables)
  
  names_real<- c("lon", "lat")
  names_pot_st<- NULL
  names_pot_end<- NULL
  sql1<- NULL
  for (i in 1:n){
    names_real<- c(names_real, paste ("filter", i, sep=""))
    names_pot_st<- c(names_pot_st, paste ("start_f", i, sep=""))
    names_pot_end<- c(names_pot_end, paste ("end_f", i, sep=""))
    sql1<- paste (sql1, paste ("real_p.filter", i, sep=""), sep=", ")   
  }
  
  names(real_p) <- names_real
  names(pot_pp) <- c(names_pot_st, names_pot_end, "groupID")
  
  filt <- names(real_p)[-c(1:2)]
  real_p$groupID <- NA
  selection_NA <- list()
  
  if (length(res) == 2) {
    for (l in 1:nrow(pot_pp)) {
      flt <-
        real_p$filter1 <= pot_pp$end_f1[l] &
        real_p$filter1 > pot_pp$start_f1[l] &
        real_p$filter2 <= pot_pp$end_f2[l] &
        real_p$filter2 > pot_pp$start_f2[l]
      real_p$groupID[flt] <- pot_pp$groupID[l]
    }
  }
  if(length(res)==3) {
    for (l in 1:nrow(pot_pp)) {
      flt <-
        real_p$filter1 <= pot_pp$end_f1[l] &
        real_p$filter1 > pot_pp$start_f1[l] &
        real_p$filter2 <= pot_pp$end_f2[l] &
        real_p$filter2 > pot_pp$start_f2[l] &
        real_p$filter3 <= pot_pp$end_f3[l] &
        real_p$filter3 > pot_pp$start_f3[l]
      real_p$groupID[flt] <- pot_pp$groupID[l]
    }
  }
  
  
  final_points <- real_p[!duplicated(real_p$groupID),]
  coord_filter<- final_points[,1:2]
  names(coord_filter)<- c("x", "y")
  
  if (plot == TRUE) {
    par (mfrow = c(1, 2), mar = c(4, 4, 0, 0.5))
    plot (
      filters[[1]],
      filters[[2]],
      pch = 19,
      col = "grey50",
      xlab = "Var 1",
      ylab = "Var 2"
    )
    points (final_points$filter1,
            final_points$filter2,
            pch = 19,
            col = "red")
    plot (coord, pch = 19, col = "grey50")
    # map(add = T)
    points (coord_filter, pch = 19, col = "red")
    par (mfrow = c(1, 1), mar = c(4, 4, 0, 0.5))
  }
  return(coord_filter)
}
