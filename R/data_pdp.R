#' Calculate data to construct partial dependence plots
#'
#' @description Calculate data to construct partial dependence plots for a given predictor
#'
#' @param model A model object of class "gam", "gbm", "glm", "graf", "ksvm", "ksvm", "maxnet”,
#' “nnet", and "randomForest" This model can be found in the first element of the list returned
#' by any function from the fit_, tune_, or esm_ function families
#' @param predictors character. Vector with a predictor name.
#' @param resolution numeric. Number of equally spaced points at which to predict continuous predictors. Default 50
#' @param resid logical. Calculate residuals based on training data. Default FALSE
#' @param training_data data.frame. Database with response (0,1) and predictor values used
#' to fit a model. Default NULL
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. When this argument is used, function will calculate partial dependence curves
#' distinguishing conditions used in training and projection conditions
#' (i.e., projection data present in projection area but not training). Default NULL
#' @param clamping logical. Perform clamping. Only for maxent models. Default FALSE
#'
#' @return A list with two tibbles "pdpdata" and "resid".
#' \itemize{
#' \item pdpdata: has data to construct partial dependence plots, the first column includes values of the selected environmental
#' variable, the second column with predicted suitability, and the third
#'  column with range type, with two values Training and Projecting, referring to suitability
#'  calculated within and outside the range of training conditions. Third column is only returned
#'  if "projection_data" argument is used
#' \item resid: has data to plot residuals. The first column includes values of the selected environmental
#'  variable and the second column with predicted suitability.
#' }
#'
#'
#' @seealso {\code{\link{data_psp}}, \code{\link{p_psp}}, \link{p_pdp}}
#'
#' @export
#' @importFrom dplyr select as_tibble
#' @importFrom gbm predict.gbm
#' @importFrom kernlab predict
#' @importFrom mgcv predict.gam
#' @importFrom stats na.omit
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
#' df <- data_pdp(
#'   model = svm_t1$model,
#'   predictors = c("aet"),
#'   resolution = 100,
#'   resid = TRUE,
#'   projection_data = somevar,
#'   training_data = abies2,
#'   clamping = FALSE
#' )
#'
#' df
#' names(df)
#' df$pdpdata
#' df$resid
#'
#' plot(df$pdpdata[1:2], type = "l")
#' points(df$resid[1:2], cex = 0.5)
#'
#' # see p_pdp to construct partial dependence plot with ggplot2
#' }
data_pdp <-
  function(model,
           predictors,
           resolution = 50,
           resid = FALSE,
           training_data = NULL,
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

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula", "ksvm", "gbm", "maxnet"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for Neural Networks, Random Forest, Support Vector Machine it is necessary to provide calibration data in 'training_data' argument"
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

    # Control average factor level
    fact <- sapply(x, is.factor)
    suit_c <- which(!fact)
    fact <- which(fact)
    suit_c <- data.frame(t(apply(x[suit_c], 2, mean)))
    # suit_c <- data.frame((x[1, ])) For residuals

    if (sum(fact) > 0) {
      for (i in 1:length(fact)) {
        ff <- sort(data.frame(unique(x[names(fact[i])]))[, 1])
        ff <- ff[as.integer(length(ff) / 2)]
        suit_c[names(fact)[i]] <- ff
      }
    }


    if (predictors %in% names(fact)) {
      rng <- sort(data.frame(unique(x[names(fact)]))[, predictors])
      suit_c <- suit_c %>% dplyr::select(-{{ predictors }})
    } else {
      if (is.null(projection_data)) {
        rng <- range(x[, predictors])
        rng <- seq(rng[1], rng[2], length.out = resolution)
      } else {
        # Range extrapolation
        rng <- terra::minmax(projection_data[[predictors]])
        rng <- seq(rng[1], rng[2], length.out = resolution)
      }
    }

    suit_c <- data.frame(rng, suit_c)
    suit_c[predictors] <- NULL
    names(suit_c)[1] <- predictors

    # Predict model
    if (class(model)[1] == "gam") {
      suit_c <-
        data.frame(suit_c[1],
                   Suitability = mgcv::predict.gam(model, newdata = suit_c, type = "response")
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = mgcv::predict.gam(model, type = "response"))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "graf") {
      suit_c <-
        data.frame(
          suit_c[1],
          Suitability = predict.graf(
            object = model,
            newdata = suit_c[names(model$peak)],
            type = "response",
            CI = NULL
          )[, 1]
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
                     Suitability = predict.graf(
                       object = model,
                       type = "response",
                       CI = NULL
                     )[, 1]
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "glm") {
      suit_c <-
        data.frame(suit_c[1],
                   Suitability = stats::predict.glm(model, newdata = suit_c, type = "response")
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = stats::predict.glm(model, type = "response"))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "gbm") {
      suit_c <-
        data.frame(suit_c[1],
                   Suitability = suppressMessages(gbm::predict.gbm(model, newdata = suit_c, type = "response"))
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = suppressMessages(gbm::predict.gbm(model, newdata = x, type = "response")))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "maxnet") {
      suit_c <-
        data.frame(suit_c[1],
                   Suitability = predict_maxnet(
                     model,
                     newdata = suit_c,
                     type = "cloglog",
                     clamp = clamping
                   )
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = predict_maxnet(
            object = model,
            newdata = data.frame(x),
            type = "cloglog",
            clamp = clamping
          ))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "nnet.formula") {
      suit_c <-
        data.frame(suit_c[1], Suitability = stats::predict(model, newdata = suit_c, type = "raw"))
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = stats::predict(model, type = "raw"))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "randomForest.formula") {
      suit_c <-
        data.frame(suit_c[1], Suitability = stats::predict(model, suit_c, type = "prob")[, 2])
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Suitability = stats::predict(model, type = "prob")[, 2])
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1],
                   Suitability = kernlab::predict(model, suit_c, type = "probabilities")[, 2]
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
                     Suitability = kernlab::predict(model, x, type = "probabilities")[, 2]
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    # Category of training and projection data
    if (!is.null(projection_data)) {
      if (!predictors %in% names(fact)) {
        result$pdpdata$Type <-
          ifelse(suit_c[, 1] >= min(x[, predictors]) &
                   suit_c[, 1] <= max(x[, predictors]),
                 "Training",
                 "Projection"
          )
      }
    }

    result <- lapply(result, function(x) if (is.data.frame(x)) dplyr::tibble(x))

    return(result)
  }
