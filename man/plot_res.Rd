% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_res.R
\name{plot_res}
\alias{plot_res}
\title{Plot different resolutions to be used in part_sblock}
\usage{
plot_res(r, res_mult)
}
\arguments{
\item{r}{SpatRaster. A raster layer, preferably a layer of environmental variables to be used}

\item{res_mult}{numeric. Maximum or minimum resolution to be tested.}
}
\value{
A plot with the original raster overlapped by a grid with the resolution used
}
\description{
This function is useful to display the maximum and minimum resolution that you want to test with the block_partition function. Note that if the resolution to be tested is very fine, the plot display may take a long time.
}
\examples{
\dontrun{
f <- system.file("external/somevar.tif", package = "flexsdm")
r <- terra::rast(f)
r <- r$CFP_1
plot_res(r, res_mult = 100)
plot_res(r, res_mult = 200)
}
}
