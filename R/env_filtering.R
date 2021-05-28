#' Function to perform environmental filtering of species occurrence records
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param variables SpatRaste. Rasters with environmental conditions
#' @param nbins integer. A number of classes used to split each environmental condition
#' @param cores integer. Number of machine cores used for processing in parallel
#'
#' @return
#' A tibble object with filtered data
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr %>% mutate select starts_with pull tibble
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats complete.cases
#' @importFrom terra extract
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Envirnomental variables
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' plot(somevar)
#'
#' # Species occurrences
#' data("spp")
#' spp
#' spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
#'
#' somevar[[1]] %>% plot()
#' points(spp1 %>% select(x, y))
#'
#' spp1$idd <- 1:nrow(spp1)
#'
#'
#' # 5 bins
#' filtered_1 <- env_filtering(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   variables = somevar,
#'   nbins = 5,
#'   cores = 1
#' )
#'
#' # 8 bins
#' filtered_2 <- env_filtering(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   variables = somevar,
#'   nbins = 8,
#'   cores = 1
#' )
#'
#' # 12 bins
#' filtered_3 <- env_filtering(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   variables = somevar,
#'   nbins = 12,
#'   cores = 1
#' )
#' # note that while higher the nbins parameter higher the number of
#' # classes to be processed (4 variables, 30 bins = 923521 classes)
#'
#' # While higher the number of bins smaller the number of records retained
#' }
env_filtering <- function(data, x, y, id, variables, nbins, cores = 1) {

  s <- . <- l <- NULL

  da <- data[c(x, y, id)]
  coord <- data[c(x, y)]

  message("Extracting values from raster ... ")
  variables <- terra::extract(variables, coord)
  variables$ID <- NULL

  filt <- stats::complete.cases(variables)
  if (sum(!filt) > 0) {
    message(sum(!filt), " records were removed because they have NAs for some variables")
    da <- da[filt, ]
    coord <- coord[filt, ]
    variables <- variables[filt, ]
  }
  rm(filt)

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
    dplyr::select(dplyr::starts_with("f")) %>%
    names()

  cl <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)

  groupID <- foreach::foreach(l = 1:nrow(real_p), .packages = c("dplyr"), .final = unlist) %dopar% {
    real_p2 <- real_p[l, ] %>% dplyr::select(dplyr::starts_with("f"))
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
  parallel::stopCluster(cl)

  final_points <- real_p[!duplicated(real_p$groupID), cnames]
  coord_filter <- da[!duplicated(real_p$groupID), c(id, x, y)]

  message("Number of filtered records: ", nrow(coord_filter))
  return(dplyr::tibble(coord_filter))
}
