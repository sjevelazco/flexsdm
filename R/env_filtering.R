#' Function to perform environmental filtering of species occurrence records
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param variables data.frame. A data.frame with environmental conditions. It is possible use two or three variables
#' @param nbins interger. A number of classes used to split each environmental condition
#' @param cores interger. Number of machine cores used for processing in parallel
env_filtering <- function(data, x, y, id, variables, nbins, cores = 1) {

  da <- data[c(x, y, id)]
  coord <- data[c(x, y)]

  message("Extracting values from raster ... ")
  variables <- raster::extract(variables, data.frame(coord))

  filt <- complete.cases(variables)
  da <- da[filt, ]
  coord <- coord[filt, ]

  variables <- variables[filt, ]
  n <- ncol(variables)
  res <- (apply(variables, 2, max) - apply(variables, 2, min)) / nbins

  classes <- list()
  for (i in 1:n) {
    ext1 <- range(variables[, i])
    ext1[1] <- ext1[1] - 1
    classes[[i]] <- seq(ext1[1], ext1[2], by = res[i])
  }
  classes <- expand.grid(classes)

  message("Number of classes in the environmental space: ", nrow(classes))
  message("Number of unfiltered records: ", nrow(da))

  ends <- NULL
  for (i in 1:n) {
    f <- classes[, i] + res[[i]]
    ends <- cbind(ends, f)
  }

  classes <- data.frame(classes, ends) %>% dplyr::mutate(groupID = c(1:nrow(classes)))
  real_p <- data.frame(coord, variables)

  names_real <- c("lon", "lat")
  names_pot_st <- NULL
  names_pot_end <- NULL
  sql1 <- NULL
  for (i in 1:n) {
    names_real <- c(names_real, paste("f", i, sep = ""))
    names_pot_st <- c(names_pot_st, paste("start_f", i, sep = ""))
    names_pot_end <- c(names_pot_end, paste("end_f", i, sep = ""))
    sql1 <- paste(sql1, paste("real_p.filter", i, sep = ""), sep = ", ")
  }

  names(real_p) <- names_real
  names(classes) <- c(names_pot_st, names_pot_end, "groupID")

  real_p$groupID <- NA

  cnames <- real_p %>%
    dplyr::select(starts_with("f")) %>%
    names()

  # pb <- progress_bar$new(total = nrow(real_p))
  # for(l in 1:nrow(real_p)) {
  #   pb$tick()
  #   real_p2 <- real_p[l,] %>% dplyr::select(starts_with('f'))
  #   flt <- list()
  #   for (ll in 1:length(cnames)) {
  #     vf <- real_p2 %>% dplyr::pull(cnames[ll])
  #     flt[[ll]] <-
  #       vf <= (classes %>% dplyr::pull(paste0('end_', cnames[ll]))) &
  #       vf > (classes %>% dplyr::pull(paste0('start_', cnames[ll])))
  #   }
  #   flt <- do.call('cbind', flt) %>% apply(., 1, all) %>% which()
  #   real_p$groupID[l] <- classes$groupID[flt]
  # }

  cl <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)

  groupID <- foreach(l = 1:nrow(real_p), .packages = c("dplyr"), .final = unlist) %dopar% {
    real_p2 <- real_p[l, ] %>% dplyr::select(starts_with("f"))
    flt <- list()
    for (ll in 1:length(cnames)) {
      vf <- real_p2 %>% dplyr::pull(cnames[ll])
      flt[[ll]] <-
        vf <= (classes %>% dplyr::pull(paste0("end_", cnames[ll]))) &
          vf > (classes %>% dplyr::pull(paste0("start_", cnames[ll])))
    }
    flt <- do.call("cbind", flt) %>%
      apply(., 1, all) %>%
      which()
    real_p$groupID[l] <- classes$groupID[flt]
    real_p$groupID[l]
  }
  real_p$groupID <- groupID
  stopCluster(cl)

  final_points <- real_p[!duplicated(real_p$groupID), cnames]
  coord_filter <- da[!duplicated(real_p$groupID), c(id, x, y)]

  message("Number of filtered records: ", nrow(coord_filter))
  return(coord_filter)
}
