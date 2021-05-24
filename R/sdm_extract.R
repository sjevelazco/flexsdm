#' Extract environmental data based on x and y coordinates
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records with x and y coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param env_layer SpatRaster. Raster with environmental variables.
#' @param variables. character. Vector with the variable names of predictor variables
#' Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers. Default NULL
#' @param filter_na logical. If filter_na = TRUE (default), the rows with NA values for any of the
#' environmental variables are removed from the returned tibble.
#'
#' @return
#'
#' A tibble that returns the original data base of presence, presence-absence, or pseudo-absence location records with additional columns
#' for the extracted environmental variables at each xy location from the SpatRast object 'env_layer'
#'
#' @importFrom dplyr tibble select all_of filter
#' @importFrom stats complete.cases
#' @importFrom terra vect crs extract
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#'
#' # Load datasets
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Extract environmental data from somevar for all species locations in spp
#' ex_spp <-
#'   sdm_extract(
#'     data = spp,
#'     x = "x",
#'     y = "y",
#'     env_layer = somevar,
#'     variables = NULL,
#'     filter_na = FALSE
#'   )
#'
#' # Extract environmental for two variables and remove rows with NAs
#' ex_spp2 <-
#'   sdm_extract(
#'     data = spp,
#'     x = "x",
#'     y = "y",
#'     env_layer = somevar,
#'     variables = c('CFP_3', 'CFP_4'),
#'     filter_na = TRUE
#'   )
#'
#' ex_spp
#' ex_spp2
#' }
sdm_extract <-
  function(data,
           x,
           y,
           env_layer,
           variables = NULL,
           filter_na = TRUE) {

    # Predictor vector when variables.=NULL
    if(is.null(variables)){
      variables <- names(env_layer)
    }

    # spatial data frame
    sp_data <-
      terra::vect(data,
                  geom = c(x, y),
                  crs = terra::crs(env_layer))

    # extract environmental data at xy locations, if filter_na = FALSE, does not remove rows with NAs
    extract_data <- dplyr::tibble(
      data,
      terra::extract(env_layer[variables],
                     sp_data,
                     cells = FALSE) %>%
        dplyr::select(dplyr::all_of(variables))
    )

    # removes rows with NAs for any environmental variable
    if (filter_na) {
      complete_vec <- stats::complete.cases(extract_data[, variables])
      if (sum(!complete_vec)>0) {
        message(sum(!complete_vec), " rows were excluded from database because NAs were found")
        extract_data <- extract_data %>% dplyr::filter(complete_vec)
      }
    }
    return(extract_data)
  }
