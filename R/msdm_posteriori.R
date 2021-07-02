#' Methods to correct overprediction of species distribution models based on occurrences and suitability patterns.
#'
#' @description These methods reduce overprediction of species distribution models based on a posteriori method (see Mendes et al 2020), i.e., the combination of the patterns of species occurrences and predicted suitability
#' @param records tibble or data.frame. A database with geographical coordinates of species presences and absences (or pseudo-absence) used to create species distribution models.
#' @param x character. Column name with longitude values.
#' @param y character. Column name with latitude values.
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
#' @param buffer character. Type o buffer width use in BMCP approach. "single" type will be used a single buffer width for all species; this value is interpreted in km (e.g. buffer=c(type="single", km=126)).
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
#' Method 'pres' (Only occurrences based restriction). This is a more restrictive variant of the OBR method. It only retains those pixels in suitability patches intercepting occurrences (k) (Mendes et al., 2020).
#'
#' Method 'lq (Lower Quantile). This method is similar to the OBR method, except by the
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
#' from Kremen et al. (2008), it is similar to the MCP except by the inclusion of a
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
#'
#' @examples
#' \dontrun{
#' require(terra)
#' data("sp_sdm") # continuous species distribution models of five species
#' data("occurrences") # presences data
#' data("absences") # absences data
#'
#' # sp_sdm is database with simple species distribution models,
#' # i.e. without any restriction method
#' plot(sp_sdm)
#'
#' # Create a temporary MSDM folder
#' tmdir <- tempdir()
#' tmdir
#' dir.create(file.path(tmdir, "MSDM"))
#' tmdir <- file.path(tmdir, "MSDM")
#' tmdir
#'
#' # The data of sp_sdm will be saved in a folder in the tmdir. This is not necessary when
#' # using your data, it is just to make this example reproducible. When you use your own data,
#' # it will be enough to have a folder with your model of the species
#' dir.create(file.path(tmdir, "original_sdm"))
#' dir_models <- file.path(tmdir, "original_sdm")
#' dir_models
#' writeRaster(sp_sdm, file.path(dir_models, names(sp_sdm)),
#'   bylayer = TRUE, format = "GTiff", overwrite = TRUE
#' )
#' # shell.exec(dir_models)
#'
#'
#' # BMCP method with a single buffer for all species----
#' MSDM_Posteriori(
#'   records = occurrences, absences = absences,
#'   x = "x", y = "y", sp = "sp", method = "BMCP", buffer = c(type = "single", km = 150),
#'   dirraster = dir_models, thr = "spec_sens",
#'   dirsave = tmdir
#' )
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#' # Categorical models corrected by BMCP methods
#' cat_bmcp <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_bmcp)
#' # Continuous models corrected by BMCP methods
#' con_bmcp <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_bmcp)
#' # shell.exec(rdir)
#'
#'
#' # OBR method----
#' MSDM_Posteriori(
#'   records = occurrences, absences = absences,
#'   x = "x", y = "y", sp = "sp", method = "OBR",
#'   dirraster = dir_models, thr = "spec_sens",
#'   dirsave = tmdir
#' )
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#'
#' # Categorical models corrected by OBR methods
#' cat_obr <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_obr)
#' # Continuous models corrected by OBR methods
#' con_obr <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_obr)
#' # shell.exec(rdir)
#' }
#'
#' @seealso \code{\link{msdm_priori}}
#' @export
msdm_posteriori <- function(records,
                            x,
                            y,
                            pr_ab,
                            method = c("OBR", "PRES", "LQ", "MCP", "BMCP"),
                            dirraster = NULL,
                            thr = "equal_sens_spec",
                            buffer = NULL) {
  if (any(is.na(c(x, y, sp)))) {
    stop("Complete 'x', 'y' or 'sp' arguments")
  }
  if (is.null(dirraster)) {
    stop("Complete 'dirraster' argument")
  }
  if (is.null(dirsave)) {
    stop("Complete 'dirsave' argument")
  }
  if (method == "BMCP" & is.null(buffer)) {
    stop("If BMCP method is used it is necessary to fill the 'buffer' argument, see the help of this function")
  }
  # if((bynarymodels==FALSE & is.null(thr))){
  #   stop("Complete 'bynarymodels' argument")
  # }
  if (any(
    thr == c(
      "kappa",
      "spec_sens",
      "no_omission",
      "prevalence",
      "equal_sens_spec",
      "sensitivty"
    )
  ) == FALSE) {
    stop(
      "'thr' argument have to be supplied with one of the next values:
      'kappa', 'spec_sens', 'no_omission',
      'prevalence', 'equal_sens_spec' or 'sensitivty'"
    )
  }


  # Create Binary folder
  foldCat <- paste(dirsave, c("BINr"), sep = "/")
  foldCon <- paste(dirsave, c("CONr"), sep = "/")
  dir.create(foldCat)
  dir.create(foldCon)

  # creation of a data.frame with presences and absences
  records <- records[, c(sp, x, y)]
  records <- records[!duplicated(records[, c(x, y)]), ]
  absences <- absences[, c(sp, x, y)]
  colnames(records) <- colnames(absences) <- c("sp", "x", "y")
  SpData <- rbind(records, absences)
  SpData$pres_abse <-
    c(rep(1, nrow(records)), rep(0, nrow(absences)))

  # Data.frame with two columns 1-names of the species
  # 2-the directory of raster of each species
  if (is.null(dirraster) == TRUE) {
    stop("Give a directory in the dirraster argument")
  } else {
    RasterList <- list.files(dirraster, pattern = ".tif$")
    sps <- gsub(".tif$", "", RasterList)
    RasterList <-
      list.files(dirraster, pattern = ".tif$", full.names = TRUE)
    RasterList <- data.frame(sps, RasterList, stringsAsFactors = FALSE)
    colnames(RasterList) <- c("sp", "RasterList")
  }

  # Vector with species names, and proving if records and raster layer
  # have the same specie names
  SpNames <- as.character(unique(records[, "sp"]))
  SpNamesR <- RasterList[, "sp"]


  if (any(!SpNames %in% SpNamesR)) {
    message(
      sum(!SpNames %in% SpNamesR),
      " species names differ between records and raster files \n"
    )
    message(
      "Next names were not found in records database: ",
      paste0(SpNamesR[!SpNames %in% SpNamesR], " ")
    )
    RasterList <- RasterList[!(RasterList[, "sp"] %in% SpNamesR[!SpNames %in% SpNamesR]), ]
    message(
      "Next names were not found in raster layers: ",
      paste0(SpNames[!SpNames %in% SpNamesR], " ")
    )
    records <- records[!(records$sp %in% SpNames[!SpNames %in% SpNamesR]), ]
    message("species names would be the same in records, absences and raster layers")
  }

  #### threshold for BMCP method
  if ((method == "BMCP" & any(buffer %in% "single"))) {
    if (all(c("type", "km") %in% names(buffer))) {
      buffer2 <- as.integer(buffer["km"]) * 1000
    } else {
      stop("Use buffer argument properly for 'single' buffer type method, e.g. buffer=c(type='single', km=100)\n")
    }
  }

  # loop to process each species
  for (s in 1:length(SpNames)) {
    cat(paste(s, "from", length(SpNames), ":", SpNames[s]), "\n")
    # Read the raster of the species
    Adeq <-
      raster::raster(RasterList[RasterList[, "sp"] == SpNames[s], "RasterList"])
    # if (is.na(crs(Adeq))) {
    #   crs(Adeq) <-
    #     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    # }

    # Extract values for one species and calculate the threshold
    singleSpData <- SpData[SpData$sp == SpNames[s], ]
    PredPoint <- raster::extract(Adeq, singleSpData[, c(x, y)])
    PredPoint <-
      data.frame(pres_abse = singleSpData[, "pres_abse"], PredPoint)
    Eval <- dismo::evaluate(
      PredPoint[PredPoint$pres_abse == 1, 2],
      PredPoint[PredPoint$pres_abse == 0, 2]
    )
    Thr <- unlist(c(dismo::threshold(Eval))[threshold])

    # MCP method----
    if (method == "MCP") {
      hull <-
        dismo::convHull(singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- dismo::predict(Adeq, hull, mask = TRUE)
      Adeq[(hull[] == 0)] <- 0
      raster::writeRaster(
        Adeq,
        paste(foldCon, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
      Mask <- Adeq > Thr
      raster::writeRaster(
        Mask,
        paste(foldCat, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
    }

    # BMCP method-----
    if (method == "BMCP") {
      hull <-
        dismo::convHull(singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- dismo::predict(Adeq, hull, mask = TRUE)
      sps <-
        singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")]
      if (raster::maxValue(hull) == 0) {
        spraster <- raster::rasterize(sps, Adeq, field = 1)
        hull[spraster[] == 1] <- 1
      }
      hull2 <- hull
      hull2[hull2[] == 0] <- NA
      hull2 <- raster::boundaries(hull2)
      hull2[hull2[] == 0] <- NA
      df <- raster::rasterToPoints(hull2)
      df <- df[df[, 3] == 1, -3]
      if (any("single" == buffer)) {
        buf <- dismo::circles(df, lonlat = TRUE, d = buffer2)
      } else if (any("species_specific" == buffer)) {
        # method based on the maximum value of the minimum distance
        dist <- flexclust::dist2(sps, sps, method = "euclidean", p = 2)
        dist[dist == 0] <- NA
        distmin <- apply(dist, 1, function(x) {
          min(x, na.rm = TRUE)
        })
        buffer2 <- sort(distmin, decreasing = T)[2]
        dist[lower.tri(dist)] <- NA
        distl <- dist == buffer2
        p1 <- which(apply(distl, 1, function(x) sum(x, na.rm = T)) == 1)
        p2 <- which(apply(distl, 2, function(x) sum(x, na.rm = T)) == 1)
        buffer2 <- distHaversine(sps[p1[1], c("x", "y")], sps[p2[1], c("x", "y")])
        # df <- (gridSample(df, aggregate(hull, fact=2) , n=1))
        buf <- circles(df, lonlat = TRUE, d = buffer2)
      }
      buf <- predict(Adeq, buf, mask = TRUE)
      buf[(hull[] == 1)] <- 1
      buf[(!is.na(Adeq[]) & is.na(buf[]))] <- 0
      Adeq[which(buf[] != 1)] <- 0

      raster::writeRaster(
        Adeq,
        paste(foldCon, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
      Mask <- (Adeq > Thr)
      raster::writeRaster(
        Mask,
        paste(foldCat, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
    }

    if (method %in% c("OBR", "LQ", "PRES")) {
      # Transform coordinate in a SpatialPoints object
      pts1 <-
        singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")]
      coordinates(pts1) <- ~ x + y
      crs(pts1) <- crs(Adeq)

      # Raster with areas equal or grater than the threshold
      AdeqBin <- Adeq >= as.numeric(Thr)
      AdeqBin[AdeqBin[] == 0] <- NA
      AdeqBin <- raster::clump(AdeqBin)
      AdeqPoints <- data.frame(raster::rasterToPoints(AdeqBin)[, 1:2])
      AdeqPoints <-
        cbind(AdeqPoints, ID = raster::extract(AdeqBin, AdeqPoints))
      # Find the patches that contain presences records
      polypoint <- as.numeric(unique(raster::extract(AdeqBin, pts1)))
      AdeqBin2 <- AdeqBin
      AdeqBin2[!AdeqBin2[] %in% polypoint] <- NA
      AdeqBin3 <- !is.na(AdeqBin2)

      # # A "SpatialPolygonsDataFrame" which each adequability patch is a feature
      # AdeqBin2 <-
      #   raster::rasterToPolygons(
      #     AdeqBin,
      #     fun = NULL,
      #     n = 8,
      #     na.rm = TRUE,
      #     digits = 12,
      #     dissolve = TRUE
      #   )
      # AdeqBin2 <- disaggregate(AdeqBin2)
      # AdeqBin2$layer <- NULL
      # # Individualize each patch with a number
      # AdeqBin2$ID <- 1:length(AdeqBin2)
      # # create a data.frame with coordinate and patch number
      # AdeqPoints <- rasterToPoints(AdeqBin)[, 1:2]
      # AdeqPoints <-
      #   cbind(AdeqPoints, ID = extract(AdeqBin2, AdeqPoints)[, 'ID'])
      # # Find the patches that contain presences records
      # polypoint <- intersect(AdeqBin2, pts1)

      # PRES methods------
      if (method == "PRES") {
        Mask <- AdeqBin3
        Mask[is.na(Mask)] <- 0
        Mask[is.na(Adeq[])] <- NA
        Mask2 <- Adeq * Mask

        raster::writeRaster(
          Mask2,
          paste(foldCon, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
        raster::writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
      } else {
        # Create a vector which contain the number (e.i. ID) of the patches
        # with presences
        filter1 <- unique(stats::na.omit(raster::values(AdeqBin2)))
        # In this step are created two data.frame one with the patches coordinates
        # that contain presences and another with patches coordinates without presences
        CoordPathP <-
          as.data.frame(AdeqPoints[AdeqPoints[, 3] %in% filter1, ])
        CoordPathNP <-
          as.data.frame(AdeqPoints[!AdeqPoints[, 3] %in% filter1, ])
        # Here is created a matrix with all combination between ID of patches
        # with and without presences

        if (ncol(CoordPathP) == 1) {
          CoordPathP <- data.frame(t(CoordPathNP))
          rownames(CoordPathP) <- NULL
        }

        if (ncol(CoordPathNP) == 1) {
          CoordPathNP <- data.frame(t(CoordPathNP))
          rownames(CoordPathNP) <- NULL
        }
        npatch1 <- unique(CoordPathP[, 3])
        npatch2 <- unique(CoordPathNP[, 3])

        DistBetweenPoly0 <- expand.grid(npatch1, npatch2)
        DistBetweenPoly0$Distance <- NA
        DistBetweenPoly0 <- as.matrix(DistBetweenPoly0)
        # Euclidean Distance between patches with and without presences
        for (i in 1:nrow(DistBetweenPoly0)) {
          comb <- (DistBetweenPoly0[i, 1:2])
          A <- CoordPathP[CoordPathP[, 3] == comb[1], 1:2]
          B <- CoordPathNP[CoordPathNP[, 3] == comb[2], 1:2]

          if (nrow(A) >= 40) {
            SEQ <- round(seq(0, nrow(A), by = (nrow(A)) / 20))
            dist <- rep(NA, length(SEQ))
            for (j in 2:length(SEQ)) {
              SEQ2 <- (SEQ[(j - 1)] + 1):SEQ[j]
              dist[j] <-
                min(flexclust::dist2(A[SEQ2, ], B, method = "euclidean", p = 2), na.rm = T)
            }
            eucdist <- min(dist[2:length(SEQ)], na.rm = T)
          } else {
            eucdist <- min(flexclust::dist2(A, B, method = "euclidean", p = 2))
          }
          DistBetweenPoly0[i, 3] <- eucdist
        }

        DistBetweenPoly0 <-
          DistBetweenPoly0[order(DistBetweenPoly0[, 2]), ]
        # Minimum Euclidean Distance between patches with and without presences
        DistBetweenPoly <-
          tapply(X = DistBetweenPoly0[, 3], DistBetweenPoly0[, 2], min)
        # # Adding value of distance patches to cells
        # AdeqBin2$Eucldist <- 0
        # AdeqBin2$Eucldist[!AdeqBin2$ID %in% filter1] <-
        #   round(DistBetweenPoly, 4)

        # OBR method------
        if (method == "OBR") {
          # method based on the maximum value of the minimum distance
          spraster <- rasterize(pts1, Adeq, field = 1)
          sps <- as(spraster, "SpatialPixels")@coords
          dist <- flexclust::dist2(sps, sps, method = "euclidean", p = 2)
          dist[dist == 0] <- NA
          distmin <- apply(dist, 1, function(x) {
            min(x, na.rm = TRUE)
          }) #
          CUT <- max(distmin)
        }
        # LQ method------
        if (method == "LQ") {
          # method based the lower quartile distance
          CUT <- c(summary(DistBetweenPoly0[, 3]))[2]
        }

        # Chosen patches
        Mask <- DistBetweenPoly0[DistBetweenPoly0[, 3] <= CUT, 2]
        Mask <-
          raster::match(AdeqBin,
                        table = c(Mask, npatch1),
                        nomatch = 0
          )
        Mask <- Mask != 0
        # Mask <- AdeqBin%in%c(Mask,npatch1)
        Mask[is.na(Adeq)] <- NA
        Mask2 <- Adeq * Mask

        # Save results as raster object
        raster::writeRaster(
          Mask2,
          paste(foldCon, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
        raster::writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], ".tif", sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
      }
    }
  }
  cat("results are in: \n", dirsave, "\n")
}
