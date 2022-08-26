#' Graphical exploration of extrapolation or suitability pattern in the environmental and geographical space
#'
#' @param training_data data.frame. Database with response (0,1) and predictor values used
#' to fit a model.
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param pr_ab character. Column name with species absence-presence, pseudo-absence-presence, or background-presence data (0,1).
#' @param extra_suit_data SpatRaster. Raster layer with extrapolation or suitability values. extra_suit_data must have the same resolution and extent than projection_data
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. projection_data must have the same resolution and extent than extra_suit_data
#' @param predictors character. Vector of predictor name(s) to calculate partial dependence plots.
#' If NULL all predictors will be used. Default NULL.
#' @param geo_space logical. If TRUE will be produced a map. Default TRUE
#' @param geo_position character. Map position regarding plot of environmental space, right, left, bottom, or upper. Default "right"
#' @param prop_points numeric. Proportion of cells from extra_suit_data and projection_data to select for plotting. default. 0.5.
#' @param maxcells integer. Maximum number of cells used to plot in the geographical space. Default 100000
#' @param alpha_p numeric. a value between 0 to 1 to control transparency of presence-absence points.
#' Lower values corresponding to more transparent colors. Default 0.5
#' @param color_p character. A vector with a color used to color presence-absence points. Default "black"
#' @param alpha_gradient numeric. a value between 0 to 1 to control transparency of projection data
#' Lower values corresponding to more transparent colors. Default 0.5
#' @param color_gradient character. A vector with colors used to color projection data. Default  c(
#' "#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")
#' @param theme ggplot2 theme. Default ggplot2::theme_classic()
#'
#' @return a plot
#' @export
#'
#' @importFrom dplyr slice_sample bind_rows
#' @importFrom ggplot2 ggplot aes_string geom_point aes scale_color_gradientn labs geom_raster scale_fill_gradientn coord_equal theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom terra is.factor as.data.frame spatSample
#' @importFrom utils combn
#'
#' @examples
#' \dontrun{
#'
#' require(dplyr)
#' require(terra)
#' require(ggplot2)
#'
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#' names(somevar) <- c("aet", "cwd", "tmx", "tmn")
#'
#' spp$species %>% unique()
#' sp <- spp %>%
#'   dplyr::filter(species == "sp2", pr_ab == 1) %>%
#'   dplyr::select(x, y, pr_ab)
#'
#' # Calibration area based on some criterion such as dispersal ability
#' ca <- calib_area(sp, x = "x", y = "y", method = c("buffer", width = 50000), crs = crs(somevar))
#'
#' plot(somevar[[1]])
#' points(sp)
#' plot(ca, add = T)
#'
#'
#' # Sampling pseudo-absences
#' set.seed(10)
#' psa <- sample_pseudoabs(
#'   data = sp,
#'   x = "x",
#'   y = "y",
#'   n = nrow(sp) * 2, # selecting number of pseudo-absence points twice number of presences
#'   method = "random",
#'   rlayer = somevar,
#'   calibarea = ca
#' )
#'
#' # Merge presences and abasences databases to get a complete calibration data
#' sp_pa <- dplyr::bind_rows(sp, psa)
#' sp_pa
#'
#' # Get environmental condition of calibration area
#' sp_pa_2 <- sdm_extract(data = sp_pa, x = "x", y = "y", env_layer = somevar)
#' sp_pa_2
#'
#' # Measure extrapolation based on calibration data (presence and pseudo-absences)
#' # using SHAPE metric
#' extr <-
#'   extra_eval(
#'     training_data = sp_pa_2, # change by training_data
#'     projection_data = somevar, # change to projection_data
#'     n_cores = 1,
#'     aggreg_factor = 1
#'   )
#' plot(extr)
#'
#' ##%######################################################%##
#' #                                                          #
#' ####            Explore extrapolation in the            ####
#' ####        environmental and geographical space        ####
#' #                                                          #
#' ##%######################################################%##
#'
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   geo_space = TRUE,
#'   prop_points = 0.05
#' )
#'
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   predictors = c("tmn", "cwd"),
#'   geo_space = TRUE,
#'   prop_points = 0.05
#' )
#'
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   predictors = c("cwd", "tmx", "aet"),
#'   geo_space = TRUE,
#'   geo_position = "left",
#'   prop_points = 0.05,
#'   color_p = "white",
#'   alpha_p = 0.5,
#'   alpha_gradient = 0.2,
#'   color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
#'   theme = ggplot2::theme_dark()
#' )
#'
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   geo_space = TRUE,
#'   prop_points = 0.05,
#'   color_p = "white",
#'   alpha_p = 0.5,
#'   alpha_gradient = 0.2,
#'   color_gradient = c("#404096", "#529DB7", "#7DB874", "#E39C37", "#D92120"),
#'   theme = ggplot2::theme_dark()
#' )
#'
#' # Explore extrapolation only in the environmental space
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   geo_space = FALSE,
#'   prop_points = 0.05,
#'   color_p = "black",
#'   color_gradient = c("#085CF8", "#65AF1E", "#F3CC1D", "#FC6A9B", "#D70500"),
#'   theme = ggplot2::theme_minimal()
#' )
#'
#'
#' ##%######################################################%##
#' #                                                          #
#' ####           With p_extra also is possible            ####
#' ####       to explore the patterns of suitability       ####
#' #                                                          #
#' ##%######################################################%##
#'
#' sp_pa_2 <- part_random(
#'   data = sp_pa_2,
#'   pr_ab = "pr_ab",
#'   method =  c(method = "kfold", folds = 5)
#' )
#'
#' rf_m1 <- fit_raf(
#'   data = sp_pa_2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmx", "tmn"),
#'   partition = ".part",
#'   thr = c("max_sorensen")
#' )
#'
#' suit <- sdm_predict(models = rf_m1, pred = somevar)
#' plot(suit$raf)
#' suit <- suit$raf
#'
#' # Pasterns of suitability in geographical and environmental space
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = suit,
#'   projection_data = somevar,
#'   geo_space = TRUE,
#'   prop_points = 0.05,
#' )
#'
#' # Pasterns of suitability selecting only presences
#' p_extra(
#'   training_data = sp_pa_2 %>%
#'     dplyr::filter(pr_ab==1),
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = suit,
#'   projection_data = somevar,
#'   geo_space = TRUE,
#'   prop_points = 0.05,
#' )
#'
#' # Pasterns of suitability in the environmental spacie only with presences
#' p_extra(
#'   training_data = sp_pa_2 %>%
#'     dplyr::filter(pr_ab==1),
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = extr,
#'   projection_data = somevar,
#'   geo_space = FALSE,
#'   prop_points = 0.05,
#' )
#' }
p_extra <- function(training_data,
                    x = "x",
                    y = "y",
                    pr_ab = "pr_ab",
                    extra_suit_data,
                    projection_data,
                    predictors = NULL,
                    geo_space = TRUE,
                    geo_position = "right",
                    prop_points = 0.2,
                    maxcells = 100000,
                    alpha_p = 0.5,
                    color_p = "black",
                    alpha_gradient = 0.5,
                    color_gradient = c(
                      "#FDE725",
                      "#B3DC2B",
                      "#6DCC57",
                      "#36B677",
                      "#1F9D87",
                      "#25818E",
                      "#30678D",
                      "#3D4988",
                      "#462777",
                      "#440154"
                    ),
                    theme = ggplot2::theme_classic()) {
  val <- Value <- NULL
  # Remove factors
  extra_suit_data <- extra_suit_data[[!terra::is.factor(extra_suit_data)]]

  if (is.null(predictors)) {
    v0 <- names(projection_data)
  } else {
    v0 <- predictors
  }
  var_comb <- utils::combn(v0, 2) %>%
    t() %>%
    data.frame()
  colnames(var_comb) <- c("x", "y")

  env_extra <- c(extra_suit_data, projection_data)
  names(env_extra)[1] <- "val"

  training_data[[pr_ab]] <- as.factor(training_data[[pr_ab]])

  p_list <- list()
  for (i in 1:nrow(var_comb)) {
    xenv <- var_comb[i, 1]
    yenv <- var_comb[i, 2]

    dfplot2 <- terra::as.data.frame(env_extra[[c(xenv, yenv, "val")]]) %>%
      unique()
    dfplot2 <- dfplot2[c(which.min(dfplot2[[3]]), which.max(dfplot2[[3]])), ]

    set.seed(10)
    dfplot <- terra::as.data.frame(env_extra[[c(xenv, yenv, "val")]]) %>%
      unique() %>%
      dplyr::slice_sample(prop = prop_points)
    dfplot <- dplyr::bind_rows(dfplot2, dfplot)

    p_list[[i]] <-
      ggplot2::ggplot(dfplot, ggplot2::aes_string(xenv, yenv)) +
      ggplot2::geom_point(ggplot2::aes(col = val), alpha = alpha_gradient) +
      ggplot2::scale_color_gradientn(colors = color_gradient) +
      ggplot2::labs(color = "Value") +
      ggplot2::geom_point(
        data = training_data,
        ggplot2::aes_string(xenv, yenv, shape = pr_ab),
        color = color_p,
        alpha = alpha_p
      )
  }
  message("Number of cell used to plot ", paste0(nrow(dfplot), " (", prop_points * 100, "%)"))


  if (geo_space) {
    rext <- terra::spatSample(extra_suit_data, maxcells, method = "regular", as.raster = TRUE) %>%
      terra::as.data.frame(rext, xy = TRUE)
    names(rext) <- c("x", "y", "Value")
    p_list[[i + 1]] <- ggplot2::ggplot(rext, ggplot2::aes(x, y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = Value)) +
      ggplot2::scale_fill_gradientn(colors = color_gradient) +
      ggplot2::geom_point(
        data = training_data,
        ggplot2::aes_string(x, y, shape = pr_ab),
        alpha = alpha_p,
        color = color_p
      )+
      ggplot2::coord_equal()
  }

  for (i in 1:length(p_list)) {
    p_list[[i]] <- p_list[[i]] + theme +
      theme(legend.title=element_blank())
  }

  if (geo_space) {
    result <- patchwork::wrap_plots(p_list[-length(p_list)]) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom") +
      theme(legend.title=element_blank())
    rmap <- (p_list[[length(p_list)]] +
               ggplot2::theme(legend.position = "none"))
    if(geo_position == "right"){
      result <- result | rmap
    } else if(geo_position == "left"){
      result <- rmap | result
    } else if(geo_position == "bottom"){
      result <- result / rmap
    } else if(geo_position == "upper"){
      result <- rmap / result
    }

  } else {
    result <- patchwork::wrap_plots(p_list) +
      patchwork::plot_layout(guides = "collect")
  }

  result
}
