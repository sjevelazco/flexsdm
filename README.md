[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
# flexsdm

### Overview 
Species distribution modeling has become a standard tool in several research areas such as ecology, 
conservation biology, biogeography, paleobiogeography, and epidemiology. Species distribution 
modeling is an area of active research in both theoretical and methodological aspects. This has 
led to defining a specific modeling process for a particular situation where the 
“click-and-run” is not a desirable option. 

One of the most exciting features of **flexsdm** is its high manipulation and parametrization 
capacity based on different functions and arguments manipulation. This enables users to define 
their own complete or partial modeling procedure specific for a modeling situation (e.g., number 
of variables, number of records, different algorithm and ensemble methods, algorithms tuning, etc.).

--- ---
### Set of functions
### 1. Pre-modeling functions 
Set of tools for preparing and manipulating modeling input data (e.g., species occurrences, 
sample pseudo-absences, define calibration area).

* `calib_area()` Delimit calibration area for constructing species distribution models
* `correct_colinvar()` Perform collinearity reduction on predictors
* `env_outliers()` Integration of outliers detection methods in the environmental space
* `part_classical()` Data partitioning for training and testing models
* `part_sblock()` Spatial block cross validation
* `plot_res()` Plot different resolutions to be used in part_sblock
* `part_senv()` Environmental and spatial cross validation
* `sample_background()` Sample background points
* `sample_pseudoabs()` Peseudo-absence allocation method
* `sdm_extract()` Extract environmental data based on x and y coordinates
* `occfilt_env()` Perform environmental filtering on species occurrences
* `occfilt_geo()` Perform geographical filtering on species occurrences


### 2. Modeling functions 
Set of function families for modeling. Several of them can be divided into `fit_*`, `tune_*`, and 
`esm_*` family functions. `fit_*` construct and validate models with default hyper-parameter 
values. `tune_*` construct and validate models searching for the best hyper-parameter values 
combination. `esm_` construct and validate Ensemble of Small Models.

#### `fit_*` functions family
* `fit_gam()` Fit and validate Generalized Additive Models
* `fit_gau()` Fit and validate Gaussian Process models
* `fit_gbm()` Fit and validate Generalized Boosted Regression models
* `fit_glm()` Fit and validate Generalized Linear Models
* `fit_max()` Fit and validate Maximum Entropy models
* `fit_net()` Fit and validate Neural Networks models
* `fit_raf()` Fit and validate Random Forest models
* `fit_svm()` Fit and validate Support Vector Machine models

#### `tune_*` functions family
* `tune_gbm()` Fit and validate Generalized Boosted Regression models with exploration of 
hyper-parameters
* `tune_max()` Fit and validate Maximum Entropy models with exploration of hyper-parameters
* `tune_net()` Fit and validate Neural Networks models with exploration of hyper-parameters
* `tune_raf()` Fit and validate Random Forest models with exploration of hyper-parameters
* `tune_svm()` Fit and validate Support Vector Machine models with exploration of hyper-parameters

#### Model ensemble
* `fit_ensemble()` Fit and validate ensemble models with different ensemble methods

#### `esm_*` functions family
* `esm_gam()` Fit and validate Generalized Additive Models with Ensemble of Small Model approach
* `esm_gau()` Fit and validate Gaussian Process models Models with Ensemble of Small Model approach
* `esm_gbm()` Fit and validate Generalized Boosted Regression models with Ensemble of Small 
Model approach
* `esm_glm()` Fit and validate Generalized Linear Models with Ensemble of Small Model approach
* `esm_max()` Fit and validate Maximum Entropy models with Ensemble of Small Model approach
* `esm_net()` Fit and validate Neural Networks models with Ensemble of Small Model approach
* `esm_svm()` Fit and validate Support Vector Machine models with Ensemble of Small Model 
approach

### 3. Post-modeling functions
Function for predicting, ensemble, and interpolate models.

* `sdm_predict()` Spatial predictions of individual and ensemble model
* `inter()` Raster interpolation between two time periods
* `extra_eval()` Measure model extrapolation
* `extra_correct()` Constraint suitability values under a given extrapolation value

--- ---
> Test the package and give us your feedback [here](https://github.com/sjevelazco/flexsdm/issues) or send an e-mail to sjevelazco@gmail.com.
