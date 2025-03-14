% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extra_eval.R
\name{extra_eval}
\alias{extra_eval}
\title{Measure model extrapolation based on Shape extrapolation metric}
\usage{
extra_eval(
  training_data,
  pr_ab,
  projection_data,
  metric = "mahalanobis",
  univar_comb = FALSE,
  n_cores = 1,
  aggreg_factor = 1
)
}
\arguments{
\item{training_data}{data.frame or tibble with environmental conditions of
presence and absence (or background points or pseudo-absences) used for constructing models}

\item{pr_ab}{character. Column name with presence and absence (or background points or
pseudo-absences) data (i.e., 1 and 0)}

\item{projection_data}{SpatRaster, data.frame or tibble with environmental condition used for projecting a model (e.g.,
a larger, encompassing region, a spatially separate region, or a different time period).
If data.frame or tibble is used function will return a tibble object.
Otherwise, as SpatRaster object.}

\item{metric}{character. Metric used to measure degree of extrapolation. Default = mahalanobis.
\itemize{
  \item mahalanobis: Degree of extrapolation is calculated based on Mahalanobis distance.
  \item euclidean: Degree of extrapolation is calculated based on Euclidean distance.
  }}

\item{univar_comb}{logical. If true, the function will add a layer or column to distinguish
between univariate (i.e., projection data outside the range of training conditions) and
combinatorial extrapolation (i.e., projection data within the range of training conditions)
using values 1 and 2, respectively. Default FALSE}

\item{n_cores}{numeric. Number of cores use for parallelization. Default 1}

\item{aggreg_factor}{positive integer. Aggregation factor expressed as number of cells in each
direction to reduce raster resolution. Use value higher than 1 would be useful when
measuring extrapolation using a raster with a high number of cells. The resolution of output will be
the same as raster object used in 'projection_data' argument. Default 1, i.e., by default, no changes
will be made to the resolution of the environmental variables.}
}
\value{
A SpatRaster or tibble object with extrapolation values measured by Shape metric. Also it
is possible estimate univariate and combinatorial extrapolation metric (see `univar_comb` argument).
}
\description{
Measure extrapolation comparing environmental data used for modeling calibration
and area for model projection. This function use the Shape metric
proposed by \href{https://doi.org/10.1111/ecog.06992}{Velazco et al., 2023}
}
\details{
This function measure model extrapolation base on the Shape metric
(\href{https://doi.org/10.1111/ecog.06992}{Velazco et al., 2023}).
Shape is a model-agnostic approach that calculates the extrapolation
degree for a given projection data point by its multivariate distance to the nearest training
data point. Such distances are relativized by a factor that reflects the dispersion of the
training data in environmental space. Distinct from other approaches (e.g.,
MESS-Multivariate Environmental Similarity Surfaces, EO-Environmental Overlap,
MOP-Mobility-Oriented Parity, EXDET-Extrapolation Detection, or AOA-Area of Applicability),
Shape incorporates an adjustable threshold to control the binary discrimination between
acceptable and unacceptable extrapolation degrees (see \code{\link{extra_truncate}}).

See this \href{https://sjevelazco.github.io/flexsdm/articles/v06_Extrapolation_example.html}{vignette at flexsdm website}
for further details about Shape metric, model truncation, and tools to explore model extrapolation.
}
\examples{
\dontrun{
require(dplyr)
require(terra)

data(spp)
f <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(f)
names(somevar) <- c("aet", "cwd", "tmx", "tmn")


spp$species \%>\% unique()
sp <- spp \%>\%
  dplyr::filter(species == "sp3", pr_ab == 1) \%>\%
  dplyr::select(x, y, pr_ab)

# Calibration area based on some criterion such as dispersal ability
ca <- calib_area(sp,
  x = "x", y = "y",
  method = c("bmcp", width = 50000),
  crs = crs(somevar)
)

plot(somevar[[1]])
points(sp)
plot(ca, add = T)


# Sampling pseudo-absences
set.seed(10)
psa <- sample_pseudoabs(
  data = sp,
  x = "x",
  y = "y",
  n = nrow(sp) * 2,
  method = "random",
  rlayer = somevar,
  calibarea = ca
)

# Merge presences and absences databases to get a complete calibration data
sp_pa <- dplyr::bind_rows(sp, psa)
sp_pa

# Get environmental condition of calibration area
sp_pa_2 <- sdm_extract(
  data = sp_pa,
  x = "x",
  y = "y",
  env_layer = somevar
)
sp_pa_2

# Measure degree of extrapolation based on Mahalanobis and
# for a projection area based on a SpatRaster object
extr <-
  extra_eval(
    training_data = sp_pa_2,
    projection_data = somevar,
    pr_ab = "pr_ab",
    n_cores = 1,
    aggreg_factor = 1,
    metric = "mahalanobis"
  )
plot(extr, main = "Extrapolation pattern")



# Let's fit, predict and truncate a model with extra_truncate
sp_pa_2 <- part_random(
  data = sp_pa_2,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

a_model <- fit_glm(
  data = sp_pa_2,
  response = "pr_ab",
  predictors = c("aet", "cwd", "tmx", "tmn"),
  partition = ".part",
  thr = c("max_sorensen")
)

predsuit <- sdm_predict(
  models = a_model,
  pred = somevar,
  thr = "max_sorensen"
)
predsuit # list with a raster with two layer
plot(predsuit[[1]])

# Truncate a model based on a given value of extrapolation
# using 'extra_truncate' function
par(mfrow = c(1, 2))
plot(extr, main = "Extrapolation")
plot(predsuit[[1]][[1]], main = "Suitability")
par(mfrow = c(1, 1))

predsuit_2 <- extra_truncate(
  suit = predsuit[[1]],
  extra = extr,
  threshold = c(50, 100, 200)
)
predsuit_2 # a list of continuous and binary models with
# different truncated at different extrapolation thresholds

plot(predsuit_2$`50`)
plot(predsuit_2$`100`)
plot(predsuit_2$`200`)


## \%######################################################\%##
####        Measure degree of extrapolation for         ####
####        projection area based on data.frame         ####
## \%######################################################\%##

extr_df <-
  extra_eval(
    training_data = sp_pa_2,
    projection_data = as.data.frame(somevar, xy = TRUE),
    pr_ab = "pr_ab",
    n_cores = 1,
    aggreg_factor = 1,
    metric = "mahalanobis"
  )
extr_df
# see 'p_extra()' to explore extrapolation or suitability pattern in the
# environmental and/or geographical space

## \%######################################################\%##
####             Explore Shape metric with              ####
####     univariate and combinatorial extrapolation     ####
## \%######################################################\%##
extr <-
  extra_eval(
    training_data = sp_pa_2,
    projection_data = somevar,
    pr_ab = "pr_ab",
    n_cores = 1,
    aggreg_factor = 1,
    metric = "mahalanobis",
    univar_comb = TRUE
  )

extr
plot(extr) # In the second layer, values equal to 1 and 2
# depict univariate and combinatorial extrapolation, respectively
}
}
\references{
\itemize{
\item Velazco, S.J.E., Brooke, M.R., De Marco Jr., P., Regan, H.M. and Franklin, J. 2023.
How far can I extrapolate my species distribution model? Exploring Shape, a novel method.
Ecography: e06992. https://doi.org/10.1111/ecog.06992
}
}
\seealso{
\code{\link{extra_truncate}}, \code{\link{p_extra}}, \code{\link{p_pdp}}, \code{\link{p_bpdp}}
}
