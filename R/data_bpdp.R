#' Calculate data to construct partial dependence surface plots
#'
#' @description Calculate data to construct Partial dependence surface plot (i.e., bivariate dependence plot) for two predictor set
#'
#' @param model A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm", "maxnet”,
#' “nnet", and "randomForest" This model can be found in the first element of the list returned
#' by any function from the fit_, tune_, or esm_ function families
#' @param predictors character. Vector with two predictor name(s) to plot. If NULL all predictors will
#' be plotted. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict continuous predictors. Default 50
#' @param training_data data.frame. Database with response (0,1) and predictor values used
#' to fit a model. Default NULL
#' @param training_boundaries character. Plot training conditions boundaries based on training
#' data (i.e., presences, presences and absences, etc).
#' If training_boundaries = "convexh", function will delimit training environmental region based on a
#' convex-hull. If training_boundaries = "rectangle", function will delimit training environmental
#' region based on four straight lines. If used any methods it is necessary provide
#' data in training_data argument.
#' If NULL all predictors will be used. Default NULL.
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. Default NULL
#' @param clamping logical. Perform clamping. Only for maxent models. Default FALSE
#'
#' @return A list with two tibbles "pdpdata" and "resid".
#' \itemize{
#' \item pspdata: has data to construct partial dependence surface plot, the first two column includes
#' values of the selected environmental variables, the third column with predicted suitability.
#' \item training_boundaries: has data to plot boundaries of training data.
#' }
#'
#'
#' @seealso \code{\link{data_pdp},  \link{p_bpdp}, \link{p_pdp}}
#'
#' @export
#' @importFrom dplyr select as_tibble
#' @importFrom gbm predict.gbm
#' @importFrom grDevices chull
#' @importFrom kernlab predict
#' @importFrom mgcv predict.gam
#' @importFrom stats na.omit predict.glm
#' @importFrom terra minmax
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
#' abies2 <- abies %>%
#'   select(x, y, pr_ab)
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
#' m <- fit_svm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmx", "tmn"),
#'   partition = ".part",
#'   thr = c("max_sens_spec")
#' )
#'
#' df <- data_bpdp(
#'   model = m$model,
#'   predictors = c("aet", "cwd"),
#'   resolution = 50,
#'   projection_data = somevar,
#'   training_boundaries = "rectangle",
#'   training_data = abies2,
#'   clamping = TRUE
#' )
#'
#' df
#' names(df)
#' df$pspdata
#' df$training_boundaries
#'
#' # see p_bpdp to construct partial dependence plot with ggplot2
#' }
data_bpdp <-
  function(model,
           predictors,
           resolution = 50,
           training_data = NULL,
           training_boundaries = NULL,
           projection_data = NULL,
           clamping = FALSE) {
    # Extract training data
    if (class(model)[1] == "gam") {
      x <- model$model[attr(model$terms, "term.labels")]
    }

    if (class(model)[1] == "graf") {
      x <- model$obsx
      x <- x[names(model$peak)]
    }

    if (class(model)[1] == "glm") {
      flt <- grepl("[I(]", attr(model$terms, "term.labels")) |
        grepl(":", attr(model$terms, "term.labels"))
      flt <- attr(model$terms, "term.labels")[!flt]
      x <- model$model[flt]
    }

    if (!is.null(training_boundaries) & is.null(training_data)) {
      stop("To extract data to delimit training boundaries it is necessary to provide training data in 'training_data' argument")
    }
    if(!is.null(training_boundaries)){
      if(!any(training_boundaries %in% c("convexh", "rectangle"))){
        stop(
          "'training_boundaries' argument could assume one of the following value: NULL, 'convexh', or 'rectangle'"
        )
      }
    }
    if (is.null(training_boundaries)) {
      training_boundaries <- 1
    }

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula", "ksvm", "gbm", "maxnet"))) {
      if (is.null(training_data)) {
        stop(
          "To estimate partial plot data for Neural Networks, Random Forest, Support Vector Machine it is necessary to provide calibration data in 'training_data' argument"
        )
      }

      if (class(model)[1] == "ksvm") {
        x <- training_data[names(attr(model@terms, "dataClasses")[-1])]
      } else if (class(model)[1] == "gbm") {
        x <- training_data[, c(model$response.name, model$var.names)]
      } else if (class(model)[1] == "maxnet") {
        x <- training_data[, names(model$samplemeans)]
      } else {
        x <- training_data[names(attr(model$terms, "dataClasses")[-1])]
      }
    }

    x <- stats::na.omit(x)
    if (training_boundaries == "convexh") {
      if(any(sapply(x[predictors], is.factor))){
        chulld <- NULL
      } else {
        chulld <- x[grDevices::chull(x[predictors]),predictors]
        chulld <- dplyr::as_tibble(chulld)
      }
    } else if (training_boundaries == "rectangle") {
      if(any(sapply(x[predictors], is.factor))){
        chulld <- NULL
      } else {
        chulld <- apply(x[predictors], 2, range)
        chulld <- expand.grid(chulld[,1], chulld[,2])
        names(chulld) <- predictors
        chulld <- dplyr::as_tibble(chulld)
      }
    } else {
      chulld <- NULL
    }


    # Control average factor level
    fact <- sapply(x, is.factor)
    suit_c <- which(!fact)
    fact <- which(fact)
    suit_c <- data.frame(t(apply(x[suit_c], 2, mean)))

    if (sum(fact) > 0) {
      for (i in 1:length(fact)) {
        ff <- sort(data.frame(unique(x[names(fact[i])]))[, 1])
        ff <- ff[as.integer(length(ff) / 2)]
        suit_c[names(fact)[i]] <- ff
      }
    }


    if (any(predictors %in% names(fact))) {
      if (is.null(projection_data)) {
        filt <- sapply(x[predictors], is.factor)
        rng1 <- range(x[predictors][!filt])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- sort(data.frame(unique(x[predictors][filt]))[,1]) # factor
      } else {
        # Range projection data
        projection_data <- projection_data[[predictors]]
        filt <- is.factor(projection_data[[predictors]])
        rng1 <- terra::minmax(projection_data[[!filt]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- as.data.frame(projection_data[[predictors]][[filt]])[,1] %>% unique()
      }
    } else {
      if (is.null(projection_data)) {
        rng1 <- range(x[, predictors[1]])
        rng2 <- range(x[, predictors[2]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- seq(rng2[1], rng2[2], length.out = resolution)
      } else {
        # Range projection data
        rng1 <- terra::minmax(projection_data[[predictors[1]]])
        rng2 <- terra::minmax(projection_data[[predictors[2]]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- seq(rng2[1], rng2[2], length.out = resolution)
      }
    }

    rng <- expand.grid(rng1, rng2)
    if(any(sapply(rng, is.factor))){
      rng <- rng[sapply(rng, is.factor)+1]
      names(rng) <- predictors[sapply(rng, is.factor)+1]
    } else {
      names(rng) <- predictors
    }
    suit_c <- suit_c %>% dplyr::select(!{{predictors}})
    suit_c <- data.frame(rng, suit_c)


    # Predict model
    if (class(model)[1] == "gam") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = mgcv::predict.gam(model, newdata = suit_c, type = "response"))
    }

    if (class(model)[1] == "graf") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = predict.graf(
                     object = model,
                     newdata = suit_c[names(model$peak)],
                     type = "response",
                     CI = NULL
                   )[, 1])
    }

    if (class(model)[1] == "glm") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = stats::predict.glm(model, newdata = suit_c, type = "response"))
    }

    if (class(model)[1] == "gbm") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = suppressMessages(gbm::predict.gbm(
                     model, newdata = suit_c, type = "response"
                   )))
    }

    if (class(model)[1] == "maxnet") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = predict_maxnet(
                     object = model,
                     newdata = suit_c,
                     type = "cloglog",
                     clamp = clamping
                   )
        )
    }

    if (class(model)[1] == "nnet.formula") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = stats::predict(model, newdata = suit_c, type = "raw"))
    }

    if (class(model)[1] == "randomForest.formula") {
      suit_c <-
        data.frame(suit_c[1:2], Suitability = stats::predict(model, suit_c, type = "prob")[, 2])
    }

    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1:2],
                   Suitability = kernlab::predict(model, suit_c, type = "probabilities")[, 2])
    }

    result <- list("pspdata" = dplyr::as_tibble(suit_c), "training_boundaries" = chulld)
    return(result)
  }
