#' Methods to correct overprediction of species distribution models based on occurrences and suitability patterns.
#'
#' @description These methods reduce overprediction of species distribution models based on a posteriori
#' methods (see Mendes et al 2020), i.e., the combination of the patterns of species occurrences
#' and predicted suitability
#'
#' @param records tibble or data.frame. A database with spatial coordinates of species presences and
#' absences (or pseudo-absence) used to create species distribution models.
#' @param x character. Column name with spatial x coordinates.
#' @param y character. Column name with spatial y coordinates.
#' @param pr_ab character. Column name with presence and absence data (i.e. 1 and 0)
#' @param method character. A character string indicating which constraint method will be used (see in details).
#' @param pres_as_patch character. For the 'lq' method, assume that cells with presences but below the threshold are patches.
#' Default FALSE.
#' @param thr character or numeric. Threshold used to get binary suitability values (i.e. 0,1), needed for
#' threshold-dependent performance metrics. Only one threshold type can be specified. It is
#' necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
#'   highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to
#'   sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   \item Also, it is possible to specify the threshold value using a numeric value (thr = 0.623)
#'   
#' Default "equal_sens_spec".
#' 
#'   }
#' 
#' @param con_thr logical. If TRUE, returns continuous suitability values for cells above the threshold, 
#' with other cells as 0. If False, returns a binary map (1 for above threshold, 0 for below).
#' Default FALSE.
#' 
#' @param buffer numeric. Buffer width use in 'bmcp' approach. The buffer width will be
#' interpreted in m if Coordinate reference system used in "crs" argument has a longitude/latitude, or map units in other
#' cases. Usage buffer=50000. Default NULL
#' @param cont_suit SpatRaster. Raster with continuous suitability predictions
#' "species_specific" type calculates the minimum pairwise-distances between all occurrences and
#' then selects the maximum distance, i.e., the value of the buffer will be the maximum distance
#' from the minimum distance. This procedure depends on the spatial pattern of each species'
#' occurrences; thus, for each species, a value of buffer width will be calculated
#' (usage buffer="species_specific").
#' @param crs character. Coordinate reference system used for calculating buffer in "bmcp" method.
#'
#' @return This function return a SpatRaster with continuous and binary prediction.
#'
#' @details
#' These function help reduce overprediction of species distribution models based on the combination
#' of the patterns of species occurrences and predicted suitability.
#' It is recommended to use these approaches only for current distribution not for models projected
#' for different time periods (past or future).
#'
#' Five methods are implemented:
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
#' These methods reduce overprediction of species distribution models already fitted
#' based on the occurrences and suitability patterns of species
#' (see 'thr' arguments)
#'
#'
#' Method 'obr' (Occurrences Based Restriction).
#' This method assumes that suitable patches intercepting species occurrences (l)
#' are more likely to be part of species distributions than suitable patches that do not
#' intercept any occurrence (k). Distance from all k patches to the closest l
#' patch is calculated, then k patches are removed that exceed a species-specific
#' distance threshold from SDMs models. This threshold (T) is calculated as the
#' maximum distance in a vector of minimum pairwise distances between occurrences.
#' Whenever a suitable pixel is within a k patch that is more than distance T from the closest l patch,
#' the suitability of the pixel is reduced to zero. It is assumed that this simple threshold
#' is a surrogate for the species-specific dispersal ability. If T is low, either the species
#' has been sampled throughout its distribution, or the species is geographically restricted,
#' justifying a narrow inclusion of k patches (Mendes et al., 2020).
#'
#' Method 'lq' (Lower Quantile). This method is similar to the 'obr' method, except by the
#' procedure to define a distance threshold to withdrawn k patches, which is the
#' lower quartile distance between k patches to the closest l patch (i.e., threshold distance is based on 
#' patches with and without presences and not in presences points as in 'obr' method) . Whenever a suitable
#' pixel is within a k patch, i.e., not within this lower quartile, the suitability of the
#' pixel is reduced to zero. This means that 75\% of k patches were withdrawn from the model (Mendes et al., 2020).
#'
#' Method 'pres' (only occurrences based restriction). This is a more restrictive variant of the 'obr' 
#' method. It only retains those pixels in suitability patches intercepting occurrences (k) (Mendes et al., 2020).
#'
#' Method 'mcp' (Minimum Convex Polygon). Compiled and adapted from
#' Kremen et al. (2008), this method excludes from SDM predictions suitable
#' pixels that do not intercept a minimum convex polygon,
#' with interior angles smaller than 180, enclosing all occurrences of a species.
#'
#' Method 'bmcp' (Buffered Minimum Convex Polygon). Compiled and adapted
#' from Kremen et al. (2008), it is similar to the 'mcp' method except for the inclusion of a
#' buffer zone surrounding minimum convex polygons. For this method a buffer width value must be provided in
#' "buffer" argument and CRS in "crs" argument.
#'
#' For further methodological and performance information of these methods see Mendes et al. (2020).
#'
#' If using one these constraining methods, cite Mendes et al (2020).
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
#' @importFrom dplyr select all_of arrange desc mutate pull filter
#' @importFrom grDevices chull
#' @importFrom methods is
#' @importFrom stats na.exclude
#' @importFrom terra rast extract vect rasterize crs buffer patches match mask unique as.polygons distance convHull
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
#' # Preparing data for modeling a species
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
#' # Fit a model
#' m_glm <- fit_glm(
#'   data = occ,
#'   response = "pr_ab",
#'   predictors = names(somevar),
#'   partition = ".part",
#'   thr = "equal_sens_spec",
#' )
#'
#'
#' # Lets predict this model
#' m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)
#' plot(m_pred[[1]])
#' m_pred[[1]] %>% plot()
#'
#' # Lets extract the raster from this list
#' m_pred <- m_pred[[1]]
#'
#' ### bmcp method
#' m_bmcp <- msdm_posteriori(
#'   records = occ,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   method = "bmcp",
#'   cont_suit = m_pred,
#'   thr = "equal_sens_spec",
#'   buffer = 30000,
#'   crs = crs(m_pred)
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
#'   cont_suit = m_pred,
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
#'   cont_suit = m_pred,
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
#'   cont_suit = m_pred,
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
#'   cont_suit = m_pred,
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
                            pres_as_patch = FALSE,
                            thr = "equal_sens_spec",
                            con_thr = FALSE,
                            buffer = NULL,
                            crs = NULL) {
  . <- thr_value <- patch <- mindis <- NULL
  method <- match.arg(method)

  if (method == "bmcp" && is.null(buffer)) {
    stop("If 'bmcp' method is used, it is necessary to fill the 'buffer' argument")
  }
  if (method == "bmcp" && is.null(crs)) {
    stop("If 'bmcp' method is used, a coordinate reference system is needed in 'crs' agument")
  }

  if (is.character(thr)) {
    valid_thr <- c("lpt", "equal_sens_spec", "max_sens_spec", "max_jaccard", "max_sorensen", "max_fpb", "sensitivity")
    if (!thr[1] %in% valid_thr) stop(paste("'thr' must be one of:", paste(valid_thr, collapse = ", ")))
  }

  #### prepare data sets
  if (!methods::is(cont_suit, "SpatRaster")) {
    cont_suit <- terra::rast(cont_suit)
  }
  cont_suit <- cont_suit[[1]]

  # database preparation
  records <- dplyr::as_tibble(records) %>%
    dplyr::select(dplyr::all_of(pr_ab), dplyr::all_of(x), dplyr::all_of(y)) %>%
    dplyr::distinct()
  colnames(records) <- c("pr_ab", "x", "y")

  # Extract values for one species and calculate the threshold
  if (!is.character(thr)) {
    thr_val <- thr[1]
  } else {
    suit_point <- terra::extract(cont_suit, records[, c("x", "y")])[[2]]
    if (thr[1] == "sensitivity") {
      thr_val <- as.numeric(thr[2])
    } else {
      eval <- sdm_eval(
        p = suit_point[records$pr_ab == 1],
        a = suit_point[records$pr_ab == 0],
        thr = thr
      )
      thr_val <- eval$thr_value[1]
    }
  }

  # Keep only presences
  pres_records <- records[records$pr_ab == 1, c("x", "y")]
  pts1 <- terra::vect(as.matrix(pres_records), crs = terra::crs(cont_suit))

  # 'mcp' method----
  if (method == "mcp") {
    hull <- terra::convHull(pts1)
    result_cont <- terra::mask(cont_suit, hull, updatevalue = 0)
  }

  # 'bmcp' method-----
  if (method == "bmcp") {
    if (identical(buffer, "species_specific")) {
      if (nrow(pts1) > 1) {
        d <- terra::distance(pts1)
        d <- as.matrix(d)
        diag(d) <- Inf
        buffer <- max(apply(d, 1, min))
      } else {
        stop("Cannot calculate 'species_specific' buffer with only one presence record")
      }
    }
    hull <- terra::convHull(pts1)
    hull <- terra::buffer(hull, width = buffer)
    result_cont <- terra::mask(cont_suit, hull, updatevalue = 0)
  }

  if (method %in% c("obr", "lq", "pres")) {
    # Raster with areas >= than the threshold
    adeq_bin <- cont_suit >= thr_val
    
    # Rasterize points to include cell without suitability but with presences points
    if(pres_as_patch && method == "lq"){
      pts1_r <- terra::rasterize(pts1, adeq_bin, background=0)
      adeq_bin <- (pts1_r + adeq_bin) > 0
    } 
    
    # Set NA to cells with 0
    adeq_bin[adeq_bin == 0] <- NA

    # Raster with patches
    patch_rast <- terra::patches(adeq_bin)
    names(patch_rast) <- "patch"

    # Find the patches that contain presences records
    patch_w_pres <- terra::extract(patch_rast, pts1)[, 2] %>%
      unique() %>%
      stats::na.exclude()
    
    # 'pres' methods------
    if (method == "pres") {
      result_cont <- terra::mask(cont_suit, patch_rast, maskvalues = patch_w_pres, inverse = TRUE, updatevalue = 0)
    } else {
      all_patches <- terra::unique(patch_rast)[, 1] %>% stats::na.exclude()
      patch_wout_pres <- setdiff(all_patches, patch_w_pres)

      poly_presence <- terra::as.polygons(patch_rast) %>% terra::subset(.$patch %in% patch_w_pres)
      poly_absence <- terra::as.polygons(patch_rast) %>% terra::subset(.$patch %in% patch_wout_pres)

      if (nrow(poly_presence) > 0 && nrow(poly_absence) > 0) {
        d_mat <- terra::distance(poly_absence, poly_presence)
        mindis_vec <- apply(d_mat, 1, min)

        # 'obr' method------
        if (method == "obr") {
          if (nrow(pts1) > 1) {
            d_pres <- terra::distance(pts1)
            d_pres <- as.matrix(d_pres)
            diag(d_pres) <- Inf
            cut <- max(apply(d_pres, 1, min))
          } else {
            cut <- 0
          }
        } else { # 'lq' method
          cut <- stats::quantile(mindis_vec, probs = 0.25)
        }

        selected_patches <- patch_wout_pres[mindis_vec <= cut]
        result_cont <- terra::mask(cont_suit, patch_rast, maskvalues = c(patch_w_pres, selected_patches), inverse = TRUE, updatevalue = 0)
      } else {
        result_cont <- terra::mask(cont_suit, patch_rast, maskvalues = patch_w_pres, inverse = TRUE, updatevalue = 0)
      }
    }
  }

  result_cont <- terra::mask(result_cont, cont_suit)

  result_bin <- result_cont >= thr_val
  result <- terra::rast(list(result_cont, result_bin))
  
  if(con_thr){
    result[[2]] <- result[[1]]*result[[2]]
  }

  names(result) <- c("msdm_cont", "msdm_bin")
  names(result)[2] <- if (is.character(thr)) thr[1] else "msdm_bin"

  return(result)
}
