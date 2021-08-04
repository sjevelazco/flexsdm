[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/sjevelazco/flexsdm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjevelazco/flexsdm/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/sjevelazco/flexsdm/branch/main/graph/badge.svg?token=UT1UB0TWSV)](https://codecov.io/gh/sjevelazco/flexsdm)

# flexsdm

### Overview 
Species distribution modeling has become a standard tool in several research areas such as ecology, conservation biology, biogeography, paleobiogeography, and epidemiology. Species distribution modeling is an area of active research in both theoretical and methodological aspects. One of the most exciting features of **flexsdm** is its high manipulation and parametrization capacity based on different functions and arguments. These attributes enable users to define a complete or partial modeling workflow specific for a modeling situation (e.g., number of variables, number of records, different algorithms, algorithms tuning, ensemble methods).

--- ---
### Structure of flexsdm
The function of **flexsdm** package are organized into three major modeling steps



<a href='https://sjevelazco.github.io/flexsdm'><img src="https://raw.githubusercontent.com/sjevelazco/flexsdm/main/man/figures/flexsdm_figure1.png" align="centre" height="550"/></a>

### 1. Pre-modeling functions 
Set tools that prepare modeling input data (e.g., species occurrences thinning, sample pseudo-absences or background points, delimitation of calibration area). 

* `calib_area()` Delimit calibration area for constructing species distribution models
* `correct_colinvar()` Collinearity reduction on predictors
* `env_outliers()` Integration of outliers detection methods in the environmental space
* `part_random()` Data partitioning for training and testing models
* `part_sblock()` Spatial block cross validation
* `part_sband()` Spatial band cross validation
* `part_senv()` Environmental cross-validation
* `plot_res()` Plot different resolutions to be used in part_sblock
* `get_block()` Transform a spatial partition layer to the same spatial properties of environmental variables
* `sample_background()` Sample background points
* `sample_pseudoabs()` Sampel pseudo-absence 
* `sdm_directory()` Create directories for saving the outputs of the flexsdm
* `sdm_extract()` Extract environmental data based on x and y coordinates
* `occfilt_env()` Perform environmental filtering on species occurrences
* `occfilt_geo()` Perform geographical filtering on species occurrences


### 2. Modeling functions 
It includes functions related to modeling construction and validation. Several of them can be grouped into `fit_*`, `tune_*`, and `esm_*` family functions. `fit_*` construct and validate models with default hyper-parameter values. `tune_*` construct and validate models searching for the best hyper-parameter values combination. `esm_` construct and validate Ensemble of Small Models.

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
Tools related to models’ geographical predictions, evaluation, and correction. 

* `sdm_predict()` Spatial predictions of individual and ensemble model
* `sdm_summarize()` Merge model performance tables
* `inter()` Raster interpolation between two time periods
* `extra_eval()` Measure model extrapolation
* `extra_correct()` Constraint suitability values under a given extrapolation value
* `msdm_priori()` Create spatial predictor variables to reduce overprediction of species distribution models
* `msdm_posteriori()` Methods to correct overprediction of species distribution models based on occurrences and suitability patterns.

-- ---
### Installation
You can install the development version of **flexsdm** from
[github](https://github.com/sjevelazco/flexsdm)

``` r
# install.packages("devtools")
devtools::install_github("sjevelazco/flexsdm")
```

--- ---
### Package website

See the package website (<https://sjevelazco.github.io/flexsdm/>) for functions explanation and vignettes.

--- ---
> Test the package and give us your feedback [here](https://github.com/sjevelazco/flexsdm/issues) or send an e-mail to sjevelazco@gmail.com.

