% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_max.R
\name{tune_max}
\alias{tune_max}
\title{Fit and validate Maximum Entropy models with exploration of hyper-parameters that optimize performance}
\usage{
tune_max(
  data,
  response,
  predictors,
  predictors_f = NULL,
  background = NULL,
  partition,
  grid = NULL,
  thr = NULL,
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 1
)
}
\arguments{
\item{data}{data.frame. Database with response (0,1) and predictors values.}

\item{response}{character. Column name with species absence-presence data (0,1).}

\item{predictors}{character. Vector with the column names of quantitative
predictor variables (i.e. continuous variables).
Usage predictors = c("aet", "cwd", "tmin")}

\item{predictors_f}{character. Vector with the column names of qualitative
predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")}

\item{background}{data.frame. Database with response variable column only containing 0 values,
and predictors variables. All
column names must be consistent with data}

\item{partition}{character. Column name with training and validation partition groups.}

\item{grid}{data.frame. A data frame object with algorithm hyper-parameters values to be tested.
It is recommended to generate this data.frame with the grid() function. Hyper-parameters needed
for tuning are 'regmult' and 'classes' (any combination of following letters l -linear-, q
 -quadratic-, h -hinge-, p -product-, and t -threshold-).}

\item{thr}{character. Threshold used to get binary suitability values (i.e. 0,1)., needed for
 threshold-dependent performance metrics. More than one threshold type can be used. It is
 necessary to provide a vector for this argument. The following threshold types are available:
\itemize{
  \item lpt: The highest threshold at which there is no omission.
  \item equal_sens_spec: Threshold at which sensitivity and specificity are equal.
  \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
  highest (aka threshold that maximizes the TSS).
  \item max_jaccard: The threshold at which the Jaccard index is the highest.
  \item max_sorensen: The threshold at which the Sorensen index is highest.
  \item max_fpb: The threshold at which # FPB (F-measure on presence-background data) is
  highest.
  \item sensitivity: Threshold based on a specified sensitivity value.
  Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to
  sensitivity value. If  a sensitivity value is not specified, a default of 0.9 will be used.
  }
If more than one threshold type is used, concatenate them, e.g., thr=c('lpt', 'max_sens_spec',
 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.}

\item{metric}{character. Performance metric used for selecting the best combination of hyper
-parameter values. One of the following metrics can be used: SORENSEN, JACCARD, FPB, TSS, KAPPA,
AUC, and BOYCE. TSS is used as default.}

\item{clamp}{logical. If TRUE, predictors and features are restricted to the range seen during
model training.}

\item{pred_type}{character. Type of response required available "link", "exponential", "cloglog"
and "logistic". Default "cloglog"}

\item{n_cores}{numeric. Number of cores use for parallelization. Default 1}
}
\value{
A list object with:
\itemize{
\item model: A "maxnet" class object from maxnet package. This object can be used for predicting.
\item predictors: A tibble with quantitative (c column names) and qualitative (f column names)
 variables use for modeling.
\item performance: Hyper-parameters values and performance metrics (see \code{\link{sdm_eval}})
for the best hyper-parameters combination.
\item hyper_performance: Performance metrics (see \code{\link{sdm_eval}}) for each combination
of the hyper-parameters.
\item data_ens: Predicted suitability for each test partition based on the best model. This
database is used in \code{\link{fit_ensemble}}
}
}
\description{
Fit and validate Maximum Entropy models with exploration of hyper-parameters that optimize performance
}
\details{
When presence-absence (or presence-pseudo-absence) data are used in data argument
in addition to background points, the function will fit models with presences and background
points and validate with presences and absences. This procedure makes maxent comparable to other
presences-absences models (e.g., random forest, support vector machine). If only presences and
background points data are used, function will fit and validate model with presences and
background data. If only presence-absences are used in data argument and without background,
function will fit model with the specified data (not recommended).
}
\examples{
\dontrun{
data("abies")
data("backg")
abies # environmental conditions of presence-absence data
backg # environmental conditions of background points

# Using k-fold partition method
# Remember that the partition method, number of folds or replications must
# be the same for presence-absence and background points datasets
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
abies2

set.seed(1)
backg <- dplyr::sample_n(backg, size = 2000, replace = FALSE)
backg2 <- part_random(
  data = backg,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
backg


gridtest <-
  expand.grid(
    regmult = seq(0.1, 3, 0.5),
    classes = c("l", "lq", "lqh")
  )

max_t1 <- tune_max(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  background = backg2,
  grid = gridtest,
  thr = "max_sens_spec",
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 2 # activate two cores to speed up this process
)

length(max_t1)
max_t1$model
max_t1$predictors
max_t1$performance
max_t1$data_ens
}
}
\seealso{
\code{\link{tune_gbm}}, \code{\link{tune_net}}, \code{\link{tune_raf}}, and \code{\link{tune_svm}}.
}
