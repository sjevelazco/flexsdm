#' Partial Dependent Plot
#'
#' @description Create partial dependence plot(s) to explore the marginal effect of
#' predictors on suitability
#'
#' @param model A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm", "maxnet”,
#' “nnet", and "randomForest" This model can be found in the first element of the list returned
#' by any function from the fit_, tune_, or esm_ function families
#' @param predictors character. Vector of predictor name(s) to calculate partial dependence plots.
#' If NULL all predictors will be used. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict suitability values for continuous
#' predictors. Default 50
#' @param resid logical. Calculate residuals based on training data. Default FALSE
#' @param training_data data.frame. Database with response (0,1) and predictor values used
#' to fit a model. Default NULL
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. When this argument is used, function will calculate partial dependence curves
#' distinguishing conditions used in training and projection conditions
#' (i.e., projection data present in projection area but not training). Default NULL
#' @param clamping logical. Perform clamping. Only for maxent models. Default FALSE
#' @param rug logical. Display training data as a rug plot on the x-axis. Note: this could be time-consuming
#'  for large databases. Default FALSE
#' @param colorl character. A vector with one or two colors used to color lines. If projection_data
#'  argument is used it is necessary to provide two colors. Default c("#462777", "#6DCC57")
#' @param colorp character. A vector with one color used to color points of residuals, Default "black"
#' @param alpha numeric. a value between 0 to 1 to control transparency of residual points.
#' Lower values corresponding to more transparent colors. Default 0.2
#' @param theme ggplot2 theme. Default ggplot2::theme_classic()
#'
#' @details This function creates partial dependent plots to explore the marginal effect of
#' predictors on suitability. If projection_data is used, function will extract the minimum and
#' maximum values found in a region or time period to which a model will be projected. If the range of projection data
#' is greater than that of the training data it will be plotted with a different color. Partial dependence curves
#' could be used to interpret a model or to explore how a model may extrapolate outside the environmental conditions
#' used to train the model.
#'
#' @seealso \code{\link{pdp_data}}, \code{\link{extra_eval}}, \code{\link{extra_truncate}}
#'
#' @importFrom ggplot2 ggplot aes_string scale_y_continuous geom_point geom_line geom_rug geom_col aes scale_color_manual geom_vline theme element_blank
#' @importFrom patchwork wrap_plots plot_layout
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
#'                       x = "x",
#'                       y = "y",
#'                       env_layer = somevar
#' )
#' abies2 <- part_random(abies2,
#'                       pr_ab = "pr_ab",
#'                       method = c(method = "kfold", folds = 5)
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
#' # Partial depence plot
#' p_pdp(model = svm_t1$model, training_data = abies2)
#' p_pdp(model = svm_t1$model, training_data = abies2, predictors = c("aet", "cwd"))
#' p_pdp(model = svm_t1$model, training_data = abies2, resolution = 5)
#' p_pdp(model = svm_t1$model, training_data = abies2, resolution = 50)
#' p_pdp(model = svm_t1$model, training_data = abies2, resid = TRUE)
#' p_pdp(model = svm_t1$model, training_data = abies2, resid = TRUE,
#'       colorl = "black", colorp = "red", alpha = 0.1)
#' p_pdp(model = svm_t1$model, training_data = abies2, resid = TRUE,
#'       colorl = "black", colorp = "red", alpha = 0.1, rug = TRUE)
#'
#' # Partial depence plot for training and projection condition found in a projection area
#' plot(somevar[[1]], main="Projection area")
#' p_pdp(model = svm_t1$model, training_data = abies2, projection_data = somevar)
#' p_pdp(model = svm_t1$model, training_data = abies2, projection_data = somevar,
#'       colorl = c("#CC00FF", "#CCFF00"))
#' p_pdp(model = svm_t1$model, training_data = abies2, projection_data = somevar,
#'       colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray")
#' p_pdp(model = svm_t1$model, training_data = abies2, projection_data = somevar,
#'       colorl = c("#CC00FF", "#CCFF00"), resid = TRUE, colorp = "gray", rug = TRUE,
#'       theme = ggplot2::theme_dark())
#' }
p_pdp <-
  function(model,
           predictors = NULL,
           resolution = 100,
           resid = FALSE,
           training_data = NULL,
           projection_data = NULL,
           clamping = FALSE,
           rug = FALSE,
           colorl = c("#462777", "#6DCC57"),
           colorp = "black",
           alpha = 0.2,
           theme = ggplot2::theme_classic()) {
    Type <- Value <- val <- NULL
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
      if (rug & is.null(training_data)) {
        stop(
          "For creating Maxent partial plot with rug it is necessary provide calibration data in 'training_data' argument"
        )
      }
      v <-
        ifelse(sapply(model$levels, is.null) == TRUE, "numeric", "factor")
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

    p <- list()

    if (is.null(projection_data)) {
      for (i in 1:length(v)) {
        crv <-
          pdp_data(
            model = model,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = NULL,
            training_data = training_data,
            clamping = clamping
          )

        if (v[i] == "numeric") {
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes_string(names(crv[[1]])[1], "Suitability")) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            {
              if (resid) ggplot2::geom_point(data = crv[[2]], color = colorp, ggplot2::aes_string(names(crv[[1]])[1], "Suitability"), alpha = alpha)
            } +
            ggplot2::geom_line(col = rev(colorl)[1], size = 0.8)

          if (rug) {
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = crv[[2]],
                ggplot2::aes_string(names(crv[[1]])[1], "Suitability"),
                sides = "b",
                alpha = 0.3
              )
          }
        } else {
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes_string(names(crv[[1]])[1], "Suitability")) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::geom_col(fill = rev(colorl)[1])
        }
      }
    } else {
      for (i in 1:length(v)) {
        crv <-
          pdp_data(
            model = model,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = projection_data[[c(names(v[i]))]],
            training_data = training_data,
            clamping = clamping
          )

        if (v[i] == "numeric") {
          rvar <- range(crv[[1]][crv[[1]]$Type == "Training", names(v[i])])

          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes_string(names(crv[[1]])[1], "Suitability")) +
            {
              if (resid) ggplot2::geom_point(data = crv[[2]], ggplot2::aes_string(names(crv[[1]])[1], "Suitability"), alpha = alpha, color = colorp)
            } +
            ggplot2::geom_line(ggplot2::aes(color = Type, group = 1), size = 0.8) +
            ggplot2::scale_color_manual(
              values = colorl,
              breaks = c("Projection", "Training"),
              name = "Range"
            ) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::geom_vline(
              xintercept = rvar,
              col = "gray70",
              linetype = 2
            )

          if (rug) {
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = crv[[2]],
                ggplot2::aes_string(names(crv[[1]])[1], "Suitability"),
                sides = "b",
                alpha = 0.5
              )
          }
        } else {
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes_string(names(crv[[1]])[1], "Suitability")) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::geom_col(fill = rev(colorl)[1])
        }
      }
    }

    # Theme
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme
    }

    # Remove y axis titles
    if (length(p) >= 4) {
      sq <- length(p) / round(sqrt(length(p)))
      sq <- seq(1, length(p), by = sq)
      sq2 <- 1:length(p)
      sq2 <- sq2[!sq2 %in% sq]
    } else if (length(p) < 4 & length(p) > 2) {
      sq2 <- 2:length(p)
    }
    if (exists("sq2")) {
      for (i in sq2) {
        p[[i]] <- p[[i]] + ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
    }

    # ncol = round(sqrt(length(p)))
    patchwork::wrap_plots(p) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  }
