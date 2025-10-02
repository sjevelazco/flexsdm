# flexsdm 1.3.x
-   `sdm_uncertainty` This new function calculates species distribution model uncertainty using a bootstrap procedure, by @sjevelazco


# flexsdm 1.3.8

-   [`fit_`](https://sjevelazco.github.io/flexsdm/index.html#fit_-functions-family), [tune\_](https://sjevelazco.github.io/flexsdm/index.html#tune_-functions-family), [esm\_](https://sjevelazco.github.io/flexsdm/index.html#esm_-functions-family) and `fit_ensemble` functions now return a performance table for each partition and replicate (performance_part), by @sjevelazco [#432](https://github.com/sjevelazco/flexsdm/pull/432) and [434](https://github.com/sjevelazco/flexsdm/pull/434)
-   `p_pdp` was improved to depict exactly training range values when projection data are used, by @sjevelazco [#429](https://github.com/sjevelazco/flexsdm/pull/429)
-   `sample_pseudoabs` a sample approach was implemented based on environmental K-means, by @sjevelazco [#410](https://github.com/sjevelazco/flexsdm/pull/410)
-   `sample_pseudoabs` K-means step was improved for the geo_env_km_const approach, by @sjevelazco [#410](https://github.com/sjevelazco/flexsdm/pull/410)
-   `p_extra` Now can handle binary maps to plot suitability or extrapolation values in the geographical or environmental space, by @sjevelazco [418](https://github.com/sjevelazco/flexsdm/pull/418)
-   `map_env_dist` This new function calculates environmental distance between presences and projection data. Only the nearest Gower distance was implemented (Domain algorithm), by @sjevelazco [419](https://github.com/sjevelazco/flexsdm/pull/419/commits/8345526cc87e50d3194030f0ef24f9202bc43a7d)
-   `fit_dom` This is a new function to fit and validate the Domain algorithm, by @sjevelazco

# flexsdm 1.3.6

-   `occfilt_geo` adapted to test different values for the three methods by @sjevelazco in [386](https://github.com/sjevelazco/flexsdm/pull/386)
-   `occfilt_env` now can filter with several bins by @sjevelazco in [388](https://github.com/sjevelazco/flexsdm/pull/388)
-   `occfilt_select` was created and a test protocol was written by @sjevelazco in [389](https://github.com/sjevelazco/flexsdm/pull/389)
-   `occfit_select` included in vignette and website by @sjevelazco in [390](https://github.com/sjevelazco/flexsdm/pull/390)
-   mean change in `sdm_predict` to handle 0 weight value for weighted mean by @sjevelazco in [391](https://github.com/sjevelazco/flexsdm/pull/391)
-   documentation of `occfilt_select` was improved by @sjevelazco in [392](https://github.com/sjevelazco/flexsdm/pull/392)
-   It is possible to tune ntree hyperparameter of random forest algorithm by @sjevelazco in [393](https://github.com/sjevelazco/flexsdm/pull/393)
-   `correct_colinvar` can be used with species points by @sjevelazco in [394](https://github.com/sjevelazco/flexsdm/pull/394)
-   Update pkg_citation.md by @mrose048 in [396](https://github.com/sjevelazco/flexsdm/pull/396)
-   `sdm_varimp` was added and tested by @mrose048 and @sjevelazco in [399](https://github.com/sjevelazco/flexsdm/pull/399)
-   min change in warning message `calib_area` by @sjevelazco in [400](https://github.com/sjevelazco/flexsdm/pull/400)

# flexsdm 1.3.4-5

-   min change in `interp` function by @sjevelazco in [362](https://github.com/sjevelazco/flexsdm/pull/362)
-   min change in `occfilt_geo` by @sjevelazco in [364](https://github.com/sjevelazco/flexsdm/pull/364)
-   Rlof was removed from dependencies by @sjevelazco in [365](https://github.com/sjevelazco/flexsdm/pull/365)
-   min change in `occfilt_geo` by @sjevelazco in [370](https://github.com/sjevelazco/flexsdm/pull/370)
-   new argument in `occfilt_geo` by @sjevelazco in [371](https://github.com/sjevelazco/flexsdm/pull/371)
-   `esm_` functions were improved by @sjevelazco in [373](https://github.com/sjevelazco/flexsdm/pull/372) & [373](https://github.com/sjevelazco/flexsdm/pull/373)
-   geographically constraint cell of env_var argument in `correct_colinvar` by @sjevelazco in [374](https://github.com/sjevelazco/flexsdm/pull/374)
-   new argument in `correct_colinvar` by @sjevelazco in [375](https://github.com/sjevelazco/flexsdm/pull/375)
-   `calib_area` was sped up by @sjevelazco in [377](https://github.com/sjevelazco/flexsdm/pull/377)
-   `correct_colivar` was fixed and improved, FA method by @sjevelazco in [383](https://github.com/sjevelazco/flexsdm/pull/383)
-   `occfilt_geo` 'defined' method was fixed by @sjevelazco in [384](https://github.com/sjevelazco/flexsdm/pull/384)

# flexsdm 1.3.4

-   it is possible to restrict the cell used to perform collinearity reduction analysis to a geographical area smaller than the full extent of environmental variables in [`correct_clinvar()`](https://sjevelazco.github.io/flexsdm/reference/correct_colinvar.html)
-   esm\_ family function was improved and debugged
-   `occfilt_geo` has a new argument "rep" to control number o repetition to filter occurrences

# flexsdm 1.3.4

-   rgeos was removed from dependencies [#356](https://github.com/sjevelazco/flexsdm/pull/356)
-   New vignette about how use different tools to explore model extrapolation and truncate models was added [#352](https://github.com/sjevelazco/flexsdm/pull/352)
-   Univariate and combinatorial extrapolation metric added to `extra_eval`.
-   Minor bugs were fixed to project PCA for other time periods [#351](https://github.com/sjevelazco/flexsdm/commit/301e241b150d75da4aa01accb3127331ca3bdcb4)
-   Best grid raster names was changed to .part in `part_sblock` and `part_sband`
-   Improvements in `correct_colinvar` to speed up function when using maxcell argument
-   Improvements in `correct_colinvar` to project PCA for other time periods

# flexsdm 1.3.3

-   Improvements in `correct_colinvar` It is now possible to sample rasters to reduce machine memory and speed up the process
-   Improvements in `sdm_predict` It is possible predict model in chunks to reduce machine memory
-   `p_extra`, `p_pdp`, and `p_bpdp` were fixed
-   New function `p_bpdp` Bivariate Partial Dependent Plot
-   New function `data_bpdp` Calculate data to construct bivariate partial dependence plots
-   Improvements in `p_dpd` Calculate data to construct partial dependence plots

# flexsdm 1.3.2

-   New function `p_extra` Graphical exploration of extrapolation or suitability pattern in the environmental and geographical space
-   New function `p_pdp` Partial Dependent Plot
-   New function `data_pdp` Calculate data to construct partial dependence plots

# flexsdm 1.3.1

-   New argument "crs" added to function `msdm_posteriori`
-   New argument "sp_name" in `sample_background` and `sample_pseudoabs`
-   raster, flexclust, ape, and sp were removed from dependencies\
-   Functions using CRS data have improved codes
-   It is possible use a numeric value to specify threshold in `msdm_posteriori`
-   `extra_eval` can use tibble or SpatRaster object in env_calib argument
-   `extra_truncate` has a new argument to define values used for model truncation and documentation was improved. \#