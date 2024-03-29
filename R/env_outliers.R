#' Integration of outliers detection methods in environmental space
#'
#' @description This function performs different methods for detecting outliers in species
#' distribution data based on the environmental conditions of occurrences. Some methods need
#' presence and absence data (e.g. Two-class Support Vector Machine and Random Forest) while other
#' only use presences (e.g. Reverse Jackknife, Box-plot, and Random Forest outliers) .
#' Outlier detection can be a useful procedure in occurrence data cleaning (Chapman 2005, Liu et al., 2018).
#'
#' @param data data.frame or tibble with presence (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data.
#' @param y character. Column name with latitude data.
#' @param id character. Column name with row id. Each row (record) must have its
#' own unique code.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param env_layer SpatRaster. Raster with environmental variables
#'
#' @details
#' This function will apply outliers detection methods to occurrence data.
#' Box-plot and Reverse Jackknife method will test outliers for each variable individually, if an
#' occurrence behaves as an outlier for at least one variable it will be highlighted as an outlier.
#' If the user uses only presence data, Support Vector Machine and Random Forest Methods will not be
#' performed. Support Vector Machine and Random Forest are performed with default
#' hyper-parameter values. In the case of a species with < 7 occurrences, the function
#' will not perform any methods (i.e. the additional columns will have 0 values); nonetheless, it will return a tibble with the additional columns with 0 and 1.
#' For further information about these methods, see Chapman (2005), Liu et al. (2018), and Velazco
#' et al. (2022).
#'
#' @return A tibble object with the same database used in 'data' argument and with seven additional columns, where 1 and 0 denote that a presence was detected or not as outliers
#' \itemize{
#'   \item .out_bxpt: outliers detected with Box-plot method
#'   \item .out_jack: outliers detected with Reverse Jackknife method
#'   \item .out_svm: outliers detected with Support Vector Machine method
#'   \item .out_rf: outliers detected with Random Forest method
#'   \item .out_rfout: outliers detected with Random Forest Outliers method
#'   \item .out_sum: frequency of a presences records was detected as outliers
#'   based on the previews methods (values between 0 and 6).
#'   }
#'
#' @references
#' \itemize{
#'   \item Chapman, A. D. (2005). Principles and methods of data cleaning: Primary Species and
#'   Species- Occurrence Data. version 1.0. Report for the Global Biodiversity Information
#'   Facility, Copenhagen. p72.  http://www.gbif.org/document/80528
#'   \item Liu, C., White, M., & Newell, G. (2018). Detecting outliers in species distribution
#'   data. Journal of Biogeography, 45(1), 164 - 176. https://doi.org/10.1111/jbi.13122
#'   \item Velazco, S.J.E.; Bedrij, N.A.; Keller, H.A.; Rojas, J.L.; Ribeiro, B.R.; De Marco, P. (2022)
#'   Quantifying the role of protected areas for safeguarding the uses of biodiversity.
#'   Biological Conservation, xx(xx) xx-xx. https://doi.org/10.1016/j.biocon.2022.109525
#'   }
#'
#' @export
#' @importFrom dplyr select mutate tibble filter pull starts_with bind_rows
#' @importFrom grDevices boxplot.stats
#' @importFrom kernlab ksvm predict
#' @importFrom randomForest randomForest outlier
#' @importFrom stats quantile
#' @importFrom terra extract
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#' require(ggplot2)
#'
#' # Environmental variables
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Species occurrences
#' data("spp")
#' spp
#' spp1 <- spp %>% dplyr::filter(species == "sp1")
#'
#' somevar[[1]] %>% plot()
#' points(spp1 %>% filter(pr_ab == 1) %>% select(x, y), col = "blue", pch = 19)
#' points(spp1 %>% filter(pr_ab == 0) %>% select(x, y), col = "red", cex = 0.5)
#'
#' spp1 <- spp1 %>% mutate(idd = 1:nrow(spp1))
#'
#' # Detect outliers
#' outs_1 <- env_outliers(
#'   data = spp1,
#'   pr_ab = "pr_ab",
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar
#' )
#'
#' # How many outliers were detected by different methods?
#' out_pa <- outs_1 %>%
#'   dplyr::select(starts_with("."), -.out_sum) %>%
#'   apply(., 2, function(x) sum(x, na.rm = T))
#' out_pa
#'
#' # How many outliers were detected by the sum of different methods?
#' outs_1 %>%
#'   dplyr::group_by(.out_sum) %>%
#'   dplyr::count()
#'
#' # Let explor where are locate records highlighted as outliers
#' outs_1 %>%
#'   dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point(aes(col = factor(.out_sum))) +
#'   facet_wrap(. ~ factor(.out_sum))
#'
#' # Detect outliers only with presences
#' outs_2 <- env_outliers(
#'   data = spp1 %>% dplyr::filter(pr_ab == 1),
#'   pr_ab = "pr_ab",
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar
#' )
#'
#' # How many outliers were detected by different methods
#' out_p <- outs_2 %>%
#'   dplyr::select(starts_with("."), -.out_sum) %>%
#'   apply(., 2, function(x) sum(x, na.rm = T))
#'
#' # How many outliers were detected by the sum of different methods?
#' outs_2 %>%
#'   dplyr::group_by(.out_sum) %>%
#'   dplyr::count()
#'
#' # Let explor where are locate records highlighted as outliers
#' outs_2 %>%
#'   dplyr::filter(pr_ab == 1, .out_sum > 0) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point(aes(col = factor(.out_sum))) +
#'   facet_wrap(. ~ factor(.out_sum))
#'
#'
#' # Comparison of function outputs when using it with
#' # presences-absences or only presences data.
#'
#' bind_rows(out_p, out_pa)
#' # Because the second case only were used presences, outliers methods
#' # based in Random Forest (.out_rf) and Support Vector Machines (.out_svm)
#' # were not performed.
#' }
env_outliers <- function(data, x, y, pr_ab, id, env_layer) {
  . <- NULL
  # Select columns and rename them
  data0 <- data
  data <- data[, c(id, x, y, pr_ab)]
  names(data) <- c("id", "x", "y", "pr_ab")

  # Convert data to tibble object
  data <- data %>% tibble()
  var <- names(env_layer)

  out_list <- list()
  occ_sp_01 <- data %>%
    dplyr::select(x, y, id, pr_ab)

  occ_sp_01 <-
    occ_sp_01 %>% dplyr::mutate(
      .out_bxpt = 0,
      .out_jack = 0,
      .out_svm = 0,
      .out_rf = 0,
      .out_rfout = 0,
      .out_sum = 0
    )

  sp_env_01 <-
    terra::extract(
      env_layer,
      terra::vect(occ_sp_01[c("x", "y")] %>%
        dplyr::rename(lon = x, lat = y))
    )[-1] %>%
    data.frame() %>%
    dplyr::tibble(id = occ_sp_01$id, pr_ab = occ_sp_01$pr_ab, .)

  # Remove NAs
  complete_vec <- stats::complete.cases(sp_env_01)
  if (sum(!complete_vec) > 0) {
    message(
      sum(!complete_vec),
      " rows were excluded from database because NAs were found"
    )
    occ_sp_01 <- occ_sp_01 %>% dplyr::filter(complete_vec)
    sp_env_01 <- sp_env_01 %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)

  sp_env_1 <- sp_env_01 %>% dplyr::filter(pr_ab == 1)
  occ_sp_01 <- occ_sp_01 %>% dplyr::select(-x, -y, -pr_ab)

  p01 <- unique(sp_env_01$pr_ab) # vector for testing presence and absence

  if (nrow(sp_env_1) > 6) {
    #### Method based on Boxplot and Reverse Jackknife ####
    l_box <- matrix(0, nrow = nrow(sp_env_1), length(var))
    l_jackk <- matrix(0, nrow = nrow(sp_env_1), length(var))
    for (ii in 1:length(var)) {
      xc <- sp_env_1 %>%
        data.frame() %>%
        dplyr::pull(var[ii])
      xr <-
        grDevices::boxplot.stats(xc, coef = 1.5)$stats %>% range()
      fe <- which((xc > xr[2] | xc < xr[1]))
      fe2 <- rev_jack(v = xc) # reverse Reverse Jackknife
      l_box[fe, ii] <- 1
      l_jackk[fe2, ii] <- 1
    }

    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_bxpt"] <-
      ifelse(rowSums(l_box) > 0, 1, 0)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_jack"] <-
      ifelse(rowSums(l_jackk) > 0, 1, 0)

    #### 	Two-class Support Vector Machine (presences absences) ####
    if (all(c(0, 1) %in% p01)) {
      sv <-
        kernlab::ksvm(
          pr_ab ~ .,
          data = sp_env_01[-1] %>% dplyr::mutate(pr_ab = factor(pr_ab)),
          type = "C-bsvc",
          kernel = "rbfdot",
          kpar = list(sigma = 0.1),
          C = 10,
          prob.model = TRUE
        )
      psv2 <-
        kernlab::predict(sv, sp_env_1, type = "probabilities")[, 2] # prediction for presences
      psv2 <- 1 - psv2 # outlierness

      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_svm"] <-
        as.integer(psv2 > stats::quantile(psv2, probs = seq(0, 1, 0.05))[20], na.rm = TRUE)
      rm(sv)

      #### Random Forest ####
      rf <-
        randomForest::randomForest(
          pr_ab ~ .,
          data = data.frame(sp_env_01[-1]) %>% dplyr::mutate(pr_ab = factor(pr_ab)),
          ntree = 2000
        )
      prd2 <- stats::predict(rf, sp_env_1[-1], "prob")[, 2]
      prd2 <- 1 - prd2 # outlierness
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rf"] <-
        as.integer(prd2 > stats::quantile(prd2, probs = seq(0, 1, 0.05))[20], na.rm = TRUE)
      rm(rf)
    }

    #### Random forest - outlier ####
    rf <-
      randomForest::randomForest(
        sp_env_1[-1] %>% dplyr::mutate(
          pr_ab =
            factor(pr_ab)
        ),
        ntree = 2000,
        proximity = TRUE
      )
    ot <- randomForest::outlier(rf)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rfout"] <-
      as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    rm(rf)

    #### Local outliers factor ####
    # if (nrow(sp_env_1) < 15) {
    #   ot <- rep(NA, nrow(sp_env_1))
    #   wii = 2
    #   while (any(is.na(ot))) {
    #     ot <- Rlof::lof(sp_env_1[-1], k = wii, cores = 1)
    #     wii <- wii + 1
    #   }
    # } else {
    #   ot <- rep(NA, nrow(sp_env_1))
    #   wii=15
    #   while(any(is.na(ot))){
    #     ot <- Rlof::lof(sp_env_1[-1], k = ifelse(wii>= nrow(sp_env_1),  nrow(sp_env_1)-1, wii), cores = 1)
    #     wii <- wii + 5
    #   }
    # }
    # occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_lof"] <-
    #   as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    # rm(ot)

    # Summary outliers
    occ_sp_01[, ".out_sum"] <-
      occ_sp_01 %>%
      dplyr::select(dplyr::starts_with(".")) %>%
      rowSums()
    out_list <- occ_sp_01
  } else {
    out_list <- occ_sp_01
  }

  out <- dplyr::bind_rows(out_list)
  cn <- "id"
  names(cn) <- id
  out <- dplyr::left_join(data0, out, by = cn)
  return(out)
}
