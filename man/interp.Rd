% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interp.R
\name{interp}
\alias{interp}
\title{Raster interpolation between two time periods}
\usage{
interp(r1, r2, y1, y2, rastername = NULL, dir_save = NULL)
}
\arguments{
\item{r1}{SpatRaster. Raster object for the initial year}

\item{r2}{SpatRaster. Raster object for the final year}

\item{y1}{numeric. Initial year}

\item{y2}{numeric. Final year}

\item{rastername}{character. Word used as prefix in raster file name. Default NULL}

\item{dir_save}{character. Directory path and name of the folder in which
the raster files will be saved. If NULL, function will return a SpatRaster object, else, it
will save raster in a given directory. Default NULL}
}
\value{
If dir_save is NULL, the function returns a SpatRaster with suitability interpolation for each year.
If dir_save is used, function outputs are saved in the directory
specified in dir_save.
}
\description{
This function interpolates values for each year between two specified years with simple interpolation using two raster objects containing e.g. habitat suitability values predicted using a species distribution model.
}
\details{
This function interpolates suitability values assuming that annual changes in
suitability are linear. This function could be useful for linking SDM output based on averaged
climate data and climate change scenarios to other models that require suitability values
disaggregated in time periods, such as population dynamics
(Keith et al., 2008; Conlisk et al., 2013; Syphard et al., 2013).
}
\examples{
\dontrun{
require(terra)
require(dplyr)

f <- system.file("external/suit_time_step.tif", package = "flexsdm")
abma <- terra::rast(f)
plot(abma)

int <- interp(
  r1 = abma[[1]],
  r2 = abma[[2]],
  y1 = 2010,
  y2 = 2020,
  rastername = "Abies",
  dir_save = NULL
)

int
}
}
\references{
\itemize{
\item Keith, D.A., Akçakaya, H.R., Thuiller, W., Midgley, G.F., Pearson, R.G., Phillips, S.J.,
 Regan, H.M., Araujo, M.B. & Rebelo, T.G. (2008) Predicting extinction risks under climate
 change: coupling stochastic population models with dynamic bioclimatic habitat models.
 Biology Letters, 4, 560-563.

\item Conlisk, E., Syphard, A.D., Franklin, J., Flint, L., Flint, A. & Regan, H.M. (2013)
 Management implications of uncertainty in assessing impacts of multiple landscape-scale
 threats to species persistence using a linked modeling approach. Global Change Biology
  3, 858-869.

\item Syphard, A.D., Regan, H.M., Franklin, J. & Swab, R. (2013) Does functional type
vulnerability to multiple threats depend on spatial context in Mediterranean-climate regions?
Diversity and Distributions, 19, 1263-1274.

}
}
