#' Extract environmental data based on x and y coordinates
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records with x and y coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param env_layer SpatRaster. Raster with environmental variables.
#' @param predictors character. Vector with the variable names of predictor variables
#' Usage predictors = c("aet", "cwd", "tmin"). If no
#' @param filter_na logical. If filter_na = TRUE (default), the rows with NA values for any of the
#' environmental variables are removed from the returned tibble.
#'
#' @return
#'
#' A tibble that returns the original data base of presence, presence-absence, or pseudo-absence location records with additional columns
#' for the extracted environmental variables at each xy location from the SpatRast object 'env_layer'
#'
#' @importFrom dplyr select all_of
#' @importFrom terra vect extract
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
#' # Extract environmental data from somevar for locations of all species in spp
#' ex_spp <-
#'   sdm_extract(
#'     data = spp %>% dplyr::filter(species == sp1),
#'     x = "x",
#'     y = "y",
#'     predictors = names(somevar),
#'     env_layer = some_var,
#'     filter_na = TRUE
#'   )
#'
#' ex_spp
#' }
sdm_extract <-
  function(data,
           x,
           y,
           predictors,
           env_layer,
           filter_na = TRUE) {
    # spatial data frame
    sp_data <-
      terra::vect(data,
                  geom = c(x, y),
                  crs = terra::crs(env_layer))

    # extract environmental data at xy locations, if filter_na = FALSE, does not remove rows with NAs
    extract_data <- data.frame(
      data,
      terra::extract(env_layer,
                     sp_data,
                     cells = TRUE) %>%
        dplyr::select(cell,
                      dplyr::all_of(predictors))
    )

    # removes rows with NAs for any environmental variable
    if (filter_na) {
      complete_vec <- stats::complete.cases(extract_data[, predictors])
    }

    extract_data <- dplyr::as_tibble(extract_data)
    return(extract_data)
  }
