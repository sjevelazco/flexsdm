#' Methods to correct overprediction of species distribution models based on occurrences and suitability patterns.
#'
#' @description These methods reduce overprediction of species distribution models based on a posteriori method (see Mendes et al 2020), i.e., the combination of the patterns of species occurrences and predicted suitability
#'
#' @param records tibble or data.frame. A database with geographical coordinates of species presences and absences (or pseudo-absence) used to create species distribution models.
#' @param x character. Column name with longitude values.
#' @param y character. Column name with latitude values.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param method character. A character string indicating which constraint method must be used create.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   }
#' In the case of use more than one threshold type it is necessary concatenate threshold types, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#' Default "equal_sens_spec".
#' @param buffer numeric. Buffer width use in 'bmcp' approach. The buffer width will be interpreted in m if raster used in cont_suit has a longitude/latitude CRS, or map units in other cases. Usage buffer=50000. Default NULL
#' @param cont_suit SpatRaster. Raster with continuous suitability predictions
#' "species_specific" type calculates the minimum pairwise-distances between all occurrences and then selects the maximum distance, i.e., the value of the buffer will be the maximum distance from the minimum distance. This procedure depends on the presence of each species' occurrences; thus, for each species, a value of buffer width will be calculated (usage buffer="species_specific").
#' @return # This function return a SpatRaster with continuous and binary prediction.
#' @details
#'
#' Abbreviation list
#'
#' \itemize{
#' \item SDM: species distribution model
#' \item l: suitability patches that intercept species occurrences
#' \item k: suitability patches that do not intercept species occurrences
#' \item T: threshold distances used to select suitability patches
#' }
#'
#'
#' These methods reduce overprediction of species distribution models already fitted
#' based on the occurrences and suitability patterns of species
#' (see 'thr' arguments)
#'
#'
#' Method 'obr' (Occurrences based restriction)-
#' This method assumes that suitable patches intercepting species occurrences (l)
#' are likely a part of species distributions than suitable patches that do not
#' intercept any occurrence (k). Distance from all patches k species occurrences to the closest l
#' patch is calculated, later it is removed k patches that trespass a species-specific
#' distance threshold from SDMs models. This threshold (T) is calculated as the
#' maximum distance in a vector of minimal pairwise distances between occurrences.
#' Whenever a suitable pixel is within a k patch distant from the closest l in less than T,
#' the suitability of the pixel was reduced to zero. We assumed that this simple threshold
#' is a surrogate of the species-specific dispersal ability. If T is low, either the species
#' has been sampled throughout its distribution, or the species is geographically restricted,
#' justifying a narrow inclusion of k patches (Mendes et al., 2020).
#'
#' Method 'pres' (Only occurrences based restriction). This is a more restrictive variant of the 'obr' method. It only retains those pixels in suitability patches intercepting occurrences (k) (Mendes et al., 2020).
#'
#' Method 'lq (Lower Quantile). This method is similar to the 'obr' method, except by the
#' procedure to define a distance threshold to withdrawn k patches, which is the
#' lower quartile distance between k patches to the closest l patch. Whenever a suitable
#' pixel is within a k patch, i.e., not within this lower quartile, the suitability of the
#' pixel is reduced to zero. This means that 75% of k patches were withdrawn from the model (Mendes et al., 2020).
#'
#' Method 'mcp' (Minimum Convex Polygon). Compiled and adapted from
#' Kremen et al. (2008), this method excludes from SDMs climate suitable
#' pixels that do not intercept a minimum convex polygon,
#' with interior angles smaller than 180, enclosing all occurrences of a species.
#'
#' Method 'bmcp' (Buffered Minimum Convex Polygon). Compiled and adapted
#' from Kremen et al. (2008), it is similar to the 'mcp' except by the inclusion of a
#' buffer zone surrounding minimum convex polygons. When used with the "single" options for buffer argument
#' function will ask for a value in km to be used as the buffer with. When used "species_specific" a buffer will be calculated for each species based on the presences occurrences patterns, assuming as buffer width
#' the maximum distance in a vector of minimal pairwise distances between occurrences.
#'
#' Further methodological and performance information of these methods see Mendes et al. (2020).
#'
#' If used one these constraining method cite Mendes et al 2020.
#'
#' @references
#' \itemize{
#' \item Mendes, P.; Velazco S.J.E.; Andrade, A.F.A.; De Marco, P. (2020) Dealing with overprediction in
#' species distribution models: how adding distance constraints can improve model accuracy,
#' Ecological Modelling, in press. https://doi.org/10.1016/j.ecolmodel.2020.109180
#' \item Kremen, C., Cameron, A., Moilanen, A., Phillips, S. J., Thomas, C. D.,
#' Beentje, H., . Zjhra, M. L. (2008). Aligning Conservation Priorities Across
#' Taxa in Madagascar with High-Resolution Planning Tools. Science, 320(5873),
#' 222-226. doi:10.1126/science.1155193
#' }
#'
#' @export
#'
#' @importFrom dplyr tibble select all_of arrange desc mutate pull filter
#' @importFrom grDevices chull
#' @importFrom sp Polygon Polygons SpatialPolygons
#' @importFrom stats na.exclude
#' @importFrom terra rast extract vect rasterize crs buffer patches match mask unique as.polygons distance
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' data("spp")
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#'
#' # It will prepared data for modeling a species
#' set.seed(10)
#' occ <- spp %>%
#'   dplyr::filter(species == "sp2") %>% # filter a species
#'   sdm_extract(
#'     data = ., x = "x", y = "y",
#'     env_layer = somevar, filter_na = TRUE
#'   ) %>% # extrac variables values
#'   part_random(.,
#'     pr_ab = "pr_ab",
#'     method = c(method = "kfold", folds = 10)
#'   ) # add columns with partition
#'
#' # Lets fit and predict a model
#' m_glm <- fit_glm(
#'   data = occ,
#'   response = "pr_ab",
#'   predictors = names(somevar),
#'   partition = ".part",
#'   thr = "equal_sens_spec",
#' )
#'
#' # Lets predict this model
#' m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)
#' plot(m_pred[[1]])
#' m_pred[[1]] %>% plot()
#'
#'
#' ### bmcp method
#' m_bmcp <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "bmcp",
#'   cont_suit = m_pred[[1]],
#'   thr = "equal_sens_spec",
#'   buffer = 30000
#' )
#'
#' plot(m_bmcp)
#'
#'
#' ### mcp method
#' m_mcp <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "mcp",
#'   cont_suit = m_pred[[1]],
#'   thr = "equal_sens_spec",
#'   buffer = NULL
#' )
#'
#' plot(m_mcp)
#'
#'
#' ### pres method
#' m_pres <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "pres",
#'   cont_suit = m_pred[[1]],
#'   thr = "equal_sens_spec",
#'   buffer = NULL
#' )
#'
#' plot(m_pres)
#'
#'
#' ### lq method
#' m_lq <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "lq",
#'   cont_suit = m_pred[[1]],
#'   thr = "equal_sens_spec",
#'   buffer = NULL
#' )
#'
#' plot(m_lq)
#'
#'
#' ### obr method
#' m_obr <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "obr",
#'   cont_suit = m_pred[[1]],
#'   thr = "equal_sens_spec",
#'   buffer = NULL
#' )
#'
#' plot(m_obr)
#' }
#'
#' @seealso \code{\link{msdm_priori}}
msdm_posteriori <- function(records,
                            x,
                            y,
                            pr_ab,
                            cont_suit,
                            method = c("obr", "pres", "lq", "mcp", "bmcp"),
                            thr = "equal_sens_spec",
                            buffer = NULL) {
  . <- thr_value <- patch <- mindis <- NULL
  if (method == "bmcp" & is.null(buffer)) {
    stop("If 'bmcp' method is used, it is necessary to fill the 'buffer' argument, see the help of this function")
  }
  if (any(
    thr == c(
      "lpt",
      "equal_sens_spec",
      "max_sens_spec",
      "max_jaccard",
      "max_sorensen",
      "max_fpb",
      "sensitivity"
    )
  ) == FALSE) {
    stop(
      "'thr' argument have to be supplied with one of the next values:
      'lpt', 'equal_sens_spec', 'max_sens_spec',
      'max_jaccard', 'max_sorensen', 'max_fpb',
      'sensitivity'"
    )
  }

  #### prepare data sets
  if (class(cont_suit) != "SpatRaster") {
    cont_suit <- terra::rast(cont_suit)
  }
  if (!any("tbl_df" %in% class(records))) {
    records <- dplyr::tibble(records)
  }

  # creation of a data.frame with presences and absences
  records <- records %>%
    dplyr::select(dplyr::all_of(pr_ab), dplyr::all_of(x), dplyr::all_of(y)) %>%
    dplyr::arrange(dplyr::all_of(dplyr::desc(pr_ab)))
  records <- records[!duplicated(records), ]
  colnames(records) <- c("pr_ab", "x", "y")

  # Extract values for one species and calculate the threshold
  suit_point <- terra::extract(cont_suit, records[, c(x, y)])[, 2]
  suit_point <-
    records %>%
    dplyr::mutate(suit_point)

  eval <-
    sdm_eval(
      p = suit_point[suit_point$pr_ab == 1, ] %>% dplyr::pull(suit_point),
      a = suit_point[suit_point$pr_ab == 0, ] %>% dplyr::pull(suit_point),
      thr = thr
    )
  thr_2 <- eval %>% dplyr::pull(thr_value)

  # 'mcp' method----
  if (method == "mcp") {
    data_pl <- data.frame(records[, c("x", "y")])
    data_pl <- data_pl[grDevices::chull(data_pl), ]
    data_pl <- sp::Polygon(data_pl)
    data_pl <- sp::Polygons(list(data_pl), ID = 1)
    data_pl <- sp::SpatialPolygons(list(data_pl))
    data_pl <- terra::vect(data_pl)
    hull <- terra::rasterize(data_pl, cont_suit)
    hull[is.na(hull)] <- 0
    result <- cont_suit * hull
    rm(hull)
    result_2 <- result >= thr_2
    result <- terra::rast(list(result, result_2))
    names(result)[2] <- thr
    return(result)
  }

  # 'bmcp' method-----
  if (method == "bmcp") {
    data_pl <- data.frame(records[, c("x", "y")])
    data_pl <- data_pl[grDevices::chull(data_pl), ]
    data_pl <- sp::Polygon(data_pl)
    data_pl <- sp::Polygons(list(data_pl), ID = 1)
    data_pl <- sp::SpatialPolygons(list(data_pl))
    data_pl <- terra::vect(data_pl)
    terra::crs(data_pl) <- terra::crs(cont_suit)
    data_pl <- terra::buffer(data_pl, width = buffer)
    hull <- terra::rasterize(data_pl, cont_suit)
    hull[is.na(hull)] <- 0
    result <- cont_suit * hull
    rm(hull)
    result_2 <- result >= thr_2
    result <- terra::rast(list(result, result_2))
    names(result)[2] <- thr
    return(result)
  }

  if (method %in% c("obr", "lq", "pres")) {
    # Transform coordinate to SpatVector object
    pts1 <- records %>%
      dplyr::filter(dplyr::all_of(pr_ab) == 1) %>%
      dplyr::select(-dplyr::all_of(pr_ab))
    pts1 <- terra::vect(as.matrix(pts1))
    terra::crs(pts1) <- terra::crs(cont_suit)

    # Raster with areas >= than the threshold
    adeq_bin <- cont_suit >= thr_2
    adeq_bin[adeq_bin == 0] <- NA

    # Raster with patches
    adeq_bin <- terra::patches(adeq_bin)
    names(adeq_bin) <- "patch"

    # Find the patches that contain presences records
    patch_w_pres <- terra::extract(adeq_bin, pts1)[, 2] %>%
      unique() %>%
      na.exclude()
    adeq_w_pres <-
      terra::match(adeq_bin, patch_w_pres) %>%
      terra::mask(adeq_bin, .)

    # 'pres' methods------
    if (method == "pres") {
      result <- cont_suit * (!is.na(adeq_w_pres))
      result_2 <- result >= thr_2
      result <- terra::rast(list(result, result_2))
      rm(result_2)
      names(result)[2] <- thr
      return(result)
    } else {
      # Create a vector which contain the number (e.i. ID) of the patches
      # with presences
      patch_w_pres <- terra::unique(adeq_w_pres)[, 1] %>% stats::na.exclude()
      patch_wout_pres <- terra::unique(adeq_bin)[, 1] %>% stats::na.exclude()
      patch_wout_pres <- patch_wout_pres[!patch_wout_pres %in% patch_w_pres]
      # In this step are created two data.frame one with the patches coordinates
      # that contain presences and another with patches coordinates without presences

      adeq_wout_np <-
        terra::match(adeq_bin, patch_wout_pres) %>%
        terra::mask(adeq_bin, .)
      poly_presence <- terra::as.polygons(adeq_w_pres)
      poly_absence <- terra::as.polygons(adeq_wout_np)

      pr_ab_poly_dist <-
        terra::distance(poly_absence, poly_presence) %>%
        data.frame() %>%
        dplyr::tibble()
      colnames(pr_ab_poly_dist) <- as.character(poly_presence$patch)
      pr_ab_poly_dist <- pr_ab_poly_dist %>%
        dplyr::mutate(patch = poly_absence$patch)
      pr_ab_poly_dist <-
        pr_ab_poly_dist %>%
        dplyr::mutate(mindis = pr_ab_poly_dist %>%
          dplyr::select(-"patch") %>% apply(., 1, min)) %>%
        dplyr::select(patch, mindis) # check if for 'lq' is used all distance or only the nearest distnace

      rm(poly_presence, poly_absence)

      # 'obr' method------
      if (method == "obr") {
        # method based on the maximum value of the minimum distance
        dist <- terra::distance(pts1) %>%
          as.matrix() %>%
          data.frame()
        dist[dist == 0] <- NA
        distmin <- apply(dist, 1, function(x) {
          min(x, na.rm = TRUE)
        }) #
        cut <- max(distmin)
      }

      # 'lq' method------
      if (method == "lq") {
        # method based the lower quartile distance
        cut <- summary(pr_ab_poly_dist$mindis)[2]
      }

      ##### Choose patches
      selected_patches <- pr_ab_poly_dist %>%
        dplyr::filter(mindis <= cut) %>%
        dplyr::pull(patch)

      filt <- terra::match(adeq_bin,
        table = c(patch_w_pres, selected_patches),
        nomatch = 0
      )
      filt <- (filt != 0)
      result <- cont_suit * filt
      result_2 <- result >= thr_2
      result <- terra::rast(list(result, result_2))
      rm(result_2, filt)
      names(result)[2] <- thr
      return(result)
    }
  }
}
