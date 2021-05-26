#' Detection of outliers in the environmental space
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param species
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param envr

env_outliers <- function(data, x, y, pr_ab, id, envr) {
  # Packages
  require(dplyr)
  require(terra)
  require(kernlab)
  require(randomForest)
  require(biogeo)
  require(Rlof)

  # Filter and replace species names
  data <- data[, c(id, x, y, pr_ab)]
  names(data) <- c("id", "species", "x", "y", "pr_ab")

  # Convert data to tibble object
  data <- data %>% tibble()
  spp <- unique(data$species)
  var <- names(envr)

  out_list <- list()
  for (i in 1:length(spp)) {
    message(i)
    occ_sp_01 <- data %>%
      dplyr::filter(species == spp[i]) %>%
      dplyr::select(species, x, y, id, pr_ab)

    occ_sp_01 <-
      occ_sp_01 %>% mutate(
        .out_bxptSUM = 0,
        .out_jackSUM = 0,
        .out_svm = 0,
        .out_rf = 0,
        .out_rfout = 0,
        .out_lof = 0,
        .out_sum = 0
      )

    sp_env_01 <-
      terra::extract(envr, data.frame(occ_sp_01[c("x", "y")])) %>%
      data.frame() %>%
      tibble(id = occ_sp_01$id, pr_ab = occ_sp_01$pr_ab, .)
    sp_env_1 <- sp_env_01 %>% dplyr::filter(pr_ab == 1)

    occ_sp_01 <- occ_sp_01 %>% dplyr::select(-x, -y, -pr_ab)

    if (nrow(sp_env_1) > 6) {
      #### Method based on Boxplot and Reverse Jackknife ####
      l_box <- matrix(0, nrow = nrow(sp_env_1), length(var))
      l_jackk <- matrix(0, nrow = nrow(sp_env_1), length(var))
      for (ii in 1:length(var)) {
        xc <- sp_env_1 %>%
          data.frame() %>%
          pull(var[ii])
        xr <- grDevices::boxplot.stats(xc, coef = 1.5)$stats %>% range()
        fe <- which((xc > xr[2] | xc < xr[1]))
        fe2 <- biogeo::rjack(xc) # reverse Reverse Jackknife
        l_box[fe, ii] <- 1
        l_jackk[fe2, ii] <- 1
      }

      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_bxptSUM"] <- ifelse(rowSums(l_box) > 0, 1, 0)
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_jackSUM"] <- ifelse(rowSums(l_jackk) > 0, 1, 0)

      #### 	Two-class Support Vector Machine (presences absences) ####
      sv <-
        kernlab::ksvm(
          pr_ab ~ .,
          data = sp_env_01[-1],
          type = "C-bsvc",
          kernel = "rbfdot",
          kpar = list(sigma = 0.1),
          C = 10,
          prob.model = T
        )
      psv2 <-
        predict(sv, sp_env_1, type = "probabilities")[, 2] # prediction for presences
      psv2 <- 1 - psv2 # outlierness

      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_svm"] <-
        as.integer(psv2 > quantile(psv2, probs = seq(0, 1, 0.05))[20])
      rm(sv)

      #### Random Forest ####
      rf <-
        randomForest::randomForest(pr_ab ~ ., data = data.frame(sp_env_01[-1]), ntree = 2000)
      prd2 <- predict(rf, sp_env_1[-1])
      prd2 <- 1 - prd2 # outlierness
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rf"] <-
        as.integer(prd2 > quantile(prd2, probs = seq(0, 1, 0.05))[20])
      rm(rf)

      #### Random forest - outlier ####
      rf <-
        randomForest::randomForest(sp_env_1[-1], ntree = 2000, proximity = TRUE)
      ot <- outlier(rf)
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_rfout"] <-
        as.integer(ot > quantile(ot, probs = seq(0, 1, 0.05))[20])
      rm(rf)

      #### Local outliers factor ####
      if (nrow(sp_env_1) < 15) {
        ot <- Rlof::lof(sp_env_1[-1], k = 5, cores = 1)
      } else {
        ot <- Rlof::lof(sp_env_1[-1], k = 15, cores = 1)
      }
      occ_sp_01[occ_sp_01$id %in% sp_env_1$id, ".out_lof"] <-
        as.integer(ot > quantile(ot, probs = seq(0, 1, 0.05))[20])
      rm(ot)

      # Summary outliers
      occ_sp_01[, ".out_sum"] <-
        occ_sp_01 %>%
        dplyr::select(starts_with(".")) %>%
        rowSums()
      out_list[[i]] <- occ_sp_01
    } else {
      out_list[[i]] <- occ_sp_01
    }
  }
  out <- dplyr::bind_rows(out_list)
  return(out)
}
