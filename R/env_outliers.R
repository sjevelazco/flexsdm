#' Several approaches for detecting outliers in the environmental space
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param id character. Column names with rows id. It is important that each row has its
#' own unique code.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param id character. Column names with rows id. It is important that each row has its
#' own unique code.
#' @param env_layer SpatRaster. Raster with environmental variables.
#'
#' @return
#' @export
#' @importFrom dplyr select mutate tibble filter pull starts_with bind_rows
#' @importFrom grDevices boxplot.stats
#' @importFrom kernlab ksvm predict
#' @importFrom randomForest randomForest outlier
#' @importFrom Rlof lof
#' @importFrom stats predict quantile
#' @importFrom terra extract
#'
#' @examples
#' \dontrun{
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' }
env_outliers <- function(data, x, y, pr_ab, id, env_layer) {
  . <- NULL
  # Select columns and rename them
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
      .out_bxptSUM = 0,
      .out_jackSUM = 0,
      .out_svm = 0,
      .out_rf = 0,
      .out_rfout = 0,
      .out_lof = 0,
      .out_sum = 0
    )

    sp_env_01 <-
    terra::extract(env_layer,
                   terra::vect(occ_sp_01[c("x", "y")] %>%
                                 dplyr::rename(lon=x, lat=y)))[-1] %>%
    data.frame() %>%
    dplyr::tibble(id = occ_sp_01$id, pr_ab = occ_sp_01$pr_ab, .)

  # Remove NAs
  complete_vec <- stats::complete.cases(sp_env_01)
  if (sum(!complete_vec) > 0) {
    message(sum(!complete_vec),
            " rows were excluded from database because NAs were found")
    occ_sp_01 <- occ_sp_01 %>% dplyr::filter(complete_vec)
    sp_env_01 <- sp_env_01 %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)

  sp_env_1 <- sp_env_01 %>% dplyr::filter(pr_ab == 1)
  occ_sp_01 <- occ_sp_01 %>% dplyr::select(-x, -y, -pr_ab)
  message('nrow sp_env_1 ', nrow(sp_env_1))
  return(sp_env_1)

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

    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_bxptSUM"] <-
      ifelse(rowSums(l_box) > 0, 1, 0)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_jackSUM"] <-
      ifelse(rowSums(l_jackk) > 0, 1, 0)

    #### 	Two-class Support Vector Machine (presences absences) ####
    if(all(c(0,1)%in%p01)){
      sv <-
        kernlab::ksvm(
          pr_ab ~ .,
          data = sp_env_01[-1] %>% dplyr::mutate(pr_ab = factor(pr_ab)),
          type = "C-bsvc",
          kernel = "rbfdot",
          kpar = list(sigma = 0.1),
          C = 10,
          prob.model = T
        )
      psv2 <-
        kernlab::predict(sv, sp_env_1, type = "probabilities")[, 2] # prediction for presences
      psv2 <- 1 - psv2 # outlierness

      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_svm"] <-
        as.integer(psv2 > stats::quantile(psv2, probs = seq(0, 1, 0.05))[20])
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
        as.integer(prd2 > stats::quantile(prd2, probs = seq(0, 1, 0.05))[20])
      rm(rf)
    }

    #### Random forest - outlier ####
    rf <-
      randomForest::randomForest(
        sp_env_1[-1] %>% dplyr::mutate(pr_ab =
                                         factor(pr_ab)),
        ntree = 2000,
        proximity = TRUE
      )
    ot <- randomForest::outlier(rf)
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rfout"] <-
      as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    rm(rf)

    #### Local outliers factor ####
    if (nrow(sp_env_1) < 15) {
      ot <- Rlof::lof(sp_env_1[-1], k = 5, cores = 1)
    } else {
      ot <- Rlof::lof(sp_env_1[-1], k = 15, cores = 1)
    }
    occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_lof"] <-
      as.integer(ot > stats::quantile(ot, probs = seq(0, 1, 0.05))[20])
    rm(ot)

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
  return(out)
}
