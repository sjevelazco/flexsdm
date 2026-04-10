##' Delimit calibration area for constructing species distribution models
#'
#' @description This function offers different methods to define the calibration area. The output could be used with other flexsdm functions like sample_backgroud, sample_pseudoabs, and sdm_predict, among others
#'
#' @param data data.frame or tibble. Database with presences
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param method character. Method used for delimiting a calibration area. Could be necessary to concatenate (c()) different objects for this argument. The following methods are implemented:
#' \itemize{
#' \item buffer: calibration area is defined by a buffer around presences. Usage method = c('buffer', width=40000). A value of the buffer width in m must be provided if CRS has a longitude/latitude, or in map units in other cases
#' \item mcp: calibration area is defined by a minimum convex polygon. Usage method = 'mcp'.
#' \item bmcp: calibration area is defined by buffered minimum convex polygon with buffer width. Usage method = c('bmcp', width=40000). A value of the buffer width in m must be provided if CRS has a longitude/latitude, or in map units in other cases
#' \item mask: calibration area is defined by selected polygons in a spatial vector object intersected by presences. Usage method = c("mask", clusters, "DN"). The second concatenated element must be a SpatVector, the third element is a character with the column name from SpatVector used for filtering polygons.
#' }
#' @param groups character. Column name indicating differentiated subsets of points. This could be used with mcp and bmcp method. Default NULL
#' @param crs character. Coordinate reference system used for transforming occurrences and outputs. If set as NULL, the result of mask method will have the same crs as the SpatVector used. Define a crs  is mandatory for buffer, mcp
#' and bmcp method.
#'
#' @return
#' A SpatVector
#' @export
#'
#' @importFrom grDevices chull
#' @importFrom methods as is
#' @importFrom terra vect buffer aggregate union crs extract is.related
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' data("spp")
#' clusters <- system.file("external/clusters.shp", package = "flexsdm")
#' clusters <- terra::vect(clusters)
#'
#'
#' single_spp <-
#'   spp %>%
#'   dplyr::filter(species == "sp1") %>%
#'   dplyr::filter(pr_ab == 1) %>%
#'   dplyr::select(-pr_ab)
#'
#'
#' plot(clusters)
#' points(single_spp[-1], col = "red")
#' crs(clusters, proj = TRUE) # coordinate reference system (CRS) used for this points database
#' # note that the unit of this CRS is in m, consequently the buffer width
#' # will be interpreted in m too
#'
#' # buffer method
#' ca_1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 40000),
#'   crs = crs(clusters)
#' )
#' plot(ca_1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mcp method
#' ca_2 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = "mcp",
#'   crs = crs(clusters)
#' )
#' plot(ca_2)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mcp method for different groups
#' single_spp <- single_spp %>% mutate(groups = ifelse(x > 150000, "a", "b"))
#'
#' plot(single_spp[, 2:3], pch = 19, col = "blue")
#' points(single_spp[single_spp$groups == "a", 2:3], col = "red", pch = 19)
#' points(single_spp[, 2:3])
#'
#' ca_2.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("mcp"),
#'   crs = crs(clusters),
#'   groups = "groups"
#' )
#' plot(ca_2.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # bmcp method
#' ca_3 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("bmcp", width = 30000),
#'   crs = crs(clusters)
#' )
#' plot(ca_3)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # bmcp method for different groups
#' ca_3.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("bmcp", width = 30000),
#'   crs = crs(clusters),
#'   groups = "groups"
#' )
#' plot(ca_3.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mask method
#' plot(clusters)
#' names(clusters)
#'
#' ca_3.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("mask", clusters, "clusters"),
#' )
#' plot(ca_3.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5, col = "red")
#' }
calib_area <- function(data, x, y, method, groups = NULL, crs = NULL) {
  . <- NULL
  m_type <- method[1]

  # Argument validation
  if (!m_type %in% c("buffer", "mcp", "bmcp", "mask")) {
    stop("argument 'method' was misused, available methods: buffer, mcp, bmcp, and mask")
  }

  if (m_type %in% c("buffer", "mcp", "bmcp") && is.null(crs)) {
    stop(paste("A coordinate reference system is needed in 'crs' argument for method:", m_type))
  }

  if (m_type %in% c("buffer", "bmcp") && !"width" %in% names(method)) {
    stop(paste0("provide width value for ", m_type, " method, e.g. method = c('", m_type, "', width = 70000)"))
  }

  if (m_type == "mask") {
    if (length(method) < 3) {
      stop("Method 'mask' requires a SpatVector and a column name, e.g. method = c('mask', vector_obj, 'col_name')")
    }
    if (!inherits(method[[2]], "SpatVector")) {
      stop("For 'mask' method, the second element of 'method' must be a SpatVector")
    }
  }

  # Data preparation
  data_cols <- c(x, y, groups)
  data <- as.data.frame(data[, data_cols])
  names(data) <- c("x", "y", "groups")[seq_along(data_cols)]

  # Execution logic
  if (m_type == "buffer") {
    data_sp <- terra::vect(data[, c("x", "y")], geom = c("x", "y"), crs = crs)
    result <- terra::buffer(data_sp, width = as.numeric(method["width"]))
    result <- terra::aggregate(result)
  }

  if (m_type %in% c("mcp", "bmcp")) {
    if (is.null(groups)) {
      data$groups <- 1
    }

    group_list <- split(data, data$groups)
    names(group_list) <- NULL
    
    result_list <- lapply(group_list, function(df) {
      # Minimum Convex Polygon
      hull_indices <- grDevices::chull(df$x, df$y)
      hull_pts <- df[hull_indices, c("x", "y")]

      # Format for terra::vect polygons (object, part, x, y, hole)
      geom_df <- data.frame(
        object = 1,
        part = 1,
        x = hull_pts$x,
        y = hull_pts$y,
        hole = 0
      )

      poly <- terra::vect(as.matrix(geom_df), type = "polygons", crs = crs)

      if (m_type == "bmcp") {
        poly <- terra::buffer(poly, width = as.numeric(method["width"]))
      }
      return(poly)
    })

    result <- if (length(result_list) > 1) {
      do.call(terra::union, result_list)
    } else {
      result_list[[1]]
    }
  }

  if (m_type == "mask") {
    polyc <- method[[2]]
    cname <- method[[3]]
    data_sp <- terra::vect(data[, c("x", "y")], geom = c("x", "y"), crs = terra::crs(polyc))

    # Select polygons in polyc that intersect any points in data_sp
    rel_matrix <- terra::is.related(polyc, data_sp, "intersects")
    if (is.matrix(rel_matrix)) {
      intersecting_rows <- rowSums(rel_matrix) > 0
    } else {
      intersecting_rows <- rel_matrix
    }
    result <- polyc[intersecting_rows, cname]
  }

  return(result)
}
