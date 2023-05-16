#' Partial Dependent Suface Plot
#'
#' @description Create partial dependence surface plot(s) to explore the bivariate marginal effect
#' of predictors on suitability
#'
#' @param model A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm", "maxnet”,
#' “nnet", and "randomForest" This model can be found in the first element of the list returned
#' by any function from the fit_, tune_, or esm_ function families
#' @param predictors character. Vector of predictor names to calculate partial dependence plots.
#' If NULL all predictors will be used. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict suitability values
#' for continuous predictors. Default 50
#' @param training_data data.frame. Database with response (0,1) and predictor values used
#' to fit a model. Default NULL
#' @param pchull logical. Plot convex-hull limit of training data. Default FALSE. If TRUE it is
#' necessary provide data in training_data argument
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. Default NULL
#' @param clamping logical. Perform clamping. Only for maxent models. Default FALSE
#' @param color_gradient character. A vector with range of colors to plot. Default c("#FDE725", "#B3DC2B",
#' "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")
#' @param color_chull character. A vector with one color used to color points of residuals, Default "white"
#' @param theme ggplot2 theme. Default ggplot2::theme_classic()
#'
#' @details This function creates partial dependent surface plots to explore the bivariate marginal
#' effect of predictors on suitability. If projection_data is used, function will extract the
#' minimum and maximum values found in a region or time period to which a model will be projected.
#' Partial dependence surface plot could be used to interpret a model or to explore how a model my
#' extrapolate outside the environmental conditions used to train the model (convex hull polygon).
#'
#' @seealso \code{\link{data_pdp}}, \code{\link{data_psp}}, \code{\link{p_pdp}},
#' \code{\link{extra_eval}}, \code{\link{extra_truncate}}
#'
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn coord_cartesian geom_polygon theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(dplyr)
#'
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar) # environmental data
#' names(somevar) <- c("aet", "cwd", "tmx", "tmn")
#' data(abies)
#'
#' # set seed
#' abies2 <- abies %>%
#'   dplyr::select(x, y, pr_ab) %>%
#'   dplyr::group_by(pr_ab) %>%
#'   dplyr::slice_sample(prop = 0.5)
#'
#' abies2 <- sdm_extract(abies2,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar
#' )
#' abies2 <- part_random(abies2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' svm_t1 <- fit_svm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmx", "tmn"),
#'   partition = ".part",
#'   thr = c("max_sens_spec")
#' )
#'
#' # Partial depence surface plot
#' p_psp(model = svm_t1$model, training_data = abies2)
#' p_psp(model = svm_t1$model, training_data = abies2, predictors = c("aet", "cwd"))
#' p_psp(model = svm_t1$model, training_data = abies2, resolution = 10)
#' p_psp(model = svm_t1$model, training_data = abies2, resolution = 70)
#' p_psp(model = svm_t1$model, training_data = abies2, pchull = TRUE)
#' p_psp(
#'   model = svm_t1$model, training_data = abies2, pchull = TRUE,
#'   color_chull = "orange",
#'   color_gradient = c("#00007F", "#007FFF", "#7FFF7F", "#FF7F00", "#7F0000")
#' )
#'
#' # Partial depence surface plot for training and projection condition
#' plot(somevar[[1]], main = "Projection area")
#' p_psp(model = svm_t1$model, training_data = abies2, projection_data = somevar, pchull = TRUE)
#'
#'
#' # PSP with categorical variables
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar) # environmental data
#' names(somevar) <- c("aet", "cwd", "tmx", "tmn")
#' cat <- system.file("external/clusters.shp", package = "flexsdm")
#' cat <- terra::vect(cat)
#' cat$clusters <- paste0("c", cat$clusters)
#' cat <- terra::rasterize(cat, somevar, field = "clusters")
#' somevar <- c(somevar, cat)
#' plot(somevar)
#'
#' # set seed
#' abies2 <- abies %>%
#'   dplyr::select(x, y, pr_ab) %>%
#'   dplyr::group_by(pr_ab) %>%
#'   dplyr::slice_sample(prop = 0.5)
#'
#' abies2 <- sdm_extract(
#'   data = abies2,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar
#' )
#' abies2 <- part_random(abies2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' svm_t1 <- fit_svm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmx", "tmn"),
#'   predictors_f = "clusters",
#'   partition = ".part",
#'   thr = c("max_sens_spec")
#' )
#'
#' p_psp(model = svm_t1$model, training_data = abies2)
#' }
p_psp <-
  function(model,
           predictors = NULL,
           resolution = 50,
           training_data = NULL,
           pchull = FALSE,
           projection_data = NULL,
           clamping = FALSE,
           color_gradient = c(
             "#000004",
             "#1B0A40",
             "#4A0C69",
             "#781B6C",
             "#A42C5F",
             "#CD4345",
             "#EC6824",
             "#FA990B",
             "#F7CF3D",
             "#FCFFA4"
           ),
           color_chull = "white",
           theme = ggplot2::theme_classic()) {
    Suitability <- Type <- Value <- val <- NULL

    if (pchull & is.null(training_data)) {
      stop(
        "For ploting partial surface plot with convex hull polygon it is necessary to provide calibration data in 'training_data' argument"
      )
    }

    if (class(model)[1] == "gam") {
      v <- attr(model$terms, "dataClasses")[-1]
    }

    if (class(model)[1] == "graf") {
      v <- sapply(model$obsx, class)
    }

    if (class(model)[1] == "glm") {
      flt <- grepl("[I(]", attr(model$terms, "term.labels")) |
        grepl(":", attr(model$terms, "term.labels"))
      flt <- attr(model$terms, "term.labels")[!flt]
      v <- attr(model$terms, "dataClasses")[flt]
    }

    if (class(model)[1] == "gbm") {
      v <- attr(model$Terms, "dataClasses")[-1]
    }

    if (class(model)[1] == "maxnet") {
      v <- ifelse(sapply(model$levels, is.null) == TRUE, "numeric", "factor")
    }

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula"))) {
      v <- attr(model$terms, "dataClasses")[-1]
    }

    if (class(model)[1] == "ksvm") {
      v <- attr(model@terms, "dataClasses")[-1]
    }

    v <- v[order(names(v))]
    if (!is.null(predictors)) {
      v <- v[names(v) %in% predictors]
    }

    var_comb <- utils::combn(names(v), 2) %>%
      t() %>%
      data.frame()

    if (any(v == "factor")) {
      filt <- (v[var_comb[, 1]] == "factor" & v[var_comb[, 2]] == "factor")
      var_comb <- var_comb[!filt, ]
      filt <- which(v[var_comb[, 1]] == "factor")
      if (length(filt) > 0) {
        var_comb[filt, ] <- var_comb[filt, 2:1]
      }
    }

    p_list <- list()

    for (i in 1:nrow(var_comb)) {
      xenv <- var_comb[i, 1]
      yenv <- var_comb[i, 2]

      crv <-
        data_psp(
          model = model,
          predictors = c(xenv, yenv),
          resolution = resolution,
          pchull = pchull,
          projection_data = projection_data,
          training_data = training_data,
          clamping = clamping
        )

      v1 <- names(crv[[1]])[1]
      v2 <- names(crv[[1]])[2]
      names(crv[[1]])[1:2] <- c("v1", "v2")

      p_list[[i]] <-
        ggplot2::ggplot(crv[[1]], ggplot2::aes(v1, v2)) +
        ggplot2::geom_raster(aes(fill = Suitability)) +
        ggplot2::scale_fill_gradientn(colours = color_gradient, limits = c(0, 1)) +
        ggplot2::coord_cartesian(expand = FALSE) +
        {
          if (!is.null(crv$pchull)) {
            names(crv[[2]]) <- c("v1", "v2")
            ggplot2::geom_polygon(
              data = crv[[2]],
              ggplot2::aes(v1, v2),
              color = color_chull,
              fill = "transparent"
            )
          }
        } +
        ggplot2::labs(x = v1, y = v2)
    }

    # Theme
    for (i in 1:length(p_list)) {
      p_list[[i]] <- p_list[[i]] + theme
    }

    if (length(p_list) == 1) {
      return(p_list[[i]])
    } else {
      result <- patchwork::wrap_plots(p_list[-length(p_list)]) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom") +
          theme(legend.title = element_blank())
      return(result)
    }
  }
