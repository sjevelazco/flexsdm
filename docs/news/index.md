# Changelog

## flexsdm 1.4.2

- `tune_max`: a grid combination whose fits all failed was assigned the
  previous combination’s model/performance (stale `mod`); now
  initialised to `NULL` and checked with
  [`is.null()`](https://rdrr.io/r/base/NULL.html)
  ([\#478](https://github.com/sjevelazco/flexsdm/issues/478))
- `fit_max`: the final refit fallback used the CV loop index in its
  check and could discard a successful fit; it now uses `is.null(mod)`
  and falls back on presences + background (`data_2`)
  ([\#478](https://github.com/sjevelazco/flexsdm/issues/478))
- `tune_svm`: removed
  [`na.omit()`](https://rspatial.github.io/terra/reference/na.omit.html)
  on partition performance so all combinations are compared across the
  same partitions
  ([\#478](https://github.com/sjevelazco/flexsdm/issues/478))
- `fit_gbm`: final model now uses `bag.fraction = 0.9`, consistent with
  the CV models; `tune_gbm` validity guard also uses 0.9
  ([\#478](https://github.com/sjevelazco/flexsdm/issues/478))
- `occfilt_geo`: moran method no longer collapses to one record (`Inf`
  distance) when no Moran’s I falls in the selection window; the closest
  value is used, and the threshold comparison is now numeric
  ([\#479](https://github.com/sjevelazco/flexsdm/issues/479))
- `sample_pseudoabs`: factor `maskval` now matches raster category IDs
  in the kmeans, env_const, env_clust and env_kmeans methods
  ([\#479](https://github.com/sjevelazco/flexsdm/issues/479))
- `get_absences`: NA in the species column no longer creates phantom NA
  rows ([\#479](https://github.com/sjevelazco/flexsdm/issues/479))
- `calib_area`: mcp/bmcp with three or more groups no longer errors
  ([\#479](https://github.com/sjevelazco/flexsdm/issues/479))
- `part_sband`: fixed uninitialised `imoran_bands_c` assignment
  ([\#479](https://github.com/sjevelazco/flexsdm/issues/479))
- CRAN review: every exported function now has a small, fast example
  that is executed on CRAN (the longer examples remain inside
  `\donttest{}`); examples of
  [`sdm_directory()`](https://sjevelazco.github.io/flexsdm/reference/sdm_directory.md)
  and
  [`homogenize_na()`](https://sjevelazco.github.io/flexsdm/reference/homogenize_na.md)
  /
  [`extra_eval()`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md)
  no longer leave files in the working directory or alter
  [`par()`](https://rdrr.io/r/graphics/par.html)
- `sdm_directory`: `main_dir` no longer defaults to a sub-folder of the
  working directory and must be provided by the user
- Functions no longer call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) with a fixed number
  (`fit_*`, `occfilt_geo`, `correct_colinvar`, `p_extra`, …). Results
  can be reproduced by calling
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before running
  them; `tune_*`, `sdm_varimp` and `sdm_uncertainty` draw the seeds used
  by parallel workers from the caller’s random-number stream
- `tune_gbm`, `tune_max`, `tune_net`, `tune_raf`, `tune_svm`,
  `sdm_varimp` and `sdm_uncertainty`: no cluster is started when
  `n_cores = 1`, the computation runs sequentially in the current
  session
- `occfilt_geo`: the `warn` option is restored when the function exits
- Gaussian process code used by
  [`fit_gau()`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md)
  and
  [`esm_gau()`](https://sjevelazco.github.io/flexsdm/reference/esm_gau.md),
  [`fit_max()`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md)
  and
  [`tune_max()`](https://sjevelazco.github.io/flexsdm/reference/tune_max.md):
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`cat()`](https://rdrr.io/r/base/cat.html) calls replaced by
  [`message()`](https://rdrr.io/r/base/message.html) or
  [`warning()`](https://rdrr.io/r/base/warning.html) (or removed)
- Check time after the last CRAN incoming check (examples took 161 s and
  tests 104 s on win-builder): the examples executed on CRAN are now
  tiny toy examples without plotting (the full examples, with their
  plots, are inside `\donttest{}`), test files that take more than a
  fraction of a second are skipped as a whole on CRAN with
  `skip_on_cran()` (they run when `NOT_CRAN=true`), and vignette `v04`
  is now precomputed (source in
  `vignettes/v04_Red_fir_example.Rmd.orig`)
- `DESCRIPTION`: reworded “ensembling” in the Description field

## flexsdm 1.4.1

- `get_absences`: Function to build a presence-absence database from
  multi-species occurrence data, by
  [@sjevelazco](https://github.com/sjevelazco)
- `occ_reproject`: Function to reproject occurrence coordinates to a new
  coordinate reference system, by
  [@sjevelazco](https://github.com/sjevelazco)
- minor warnings and package formality were fixed, by
  [@sjevelazco](https://github.com/sjevelazco)
  [467](https://github.com/sjevelazco/flexsdm/pull/467)
- `correct_colinvar`: Unicode box-drawing characters in the
  documentation’s example directory tree were replaced with plain ASCII,
  fixing a LaTeX error that caused the PDF manual build to fail on CRAN,
  by [@sjevelazco](https://github.com/sjevelazco)
- Vignette `v07_Complete_workflow`: two dataset links were missing the
  `https://` scheme (`www.doi.org/...`), which CRAN flagged as invalid
  URIs; fixed to `https://doi.org/...`, by
  [@sjevelazco](https://github.com/sjevelazco)
- `tune_max`: fixed a bug where the default hyperparameter grid (used
  when `grid = NULL`) tested a `regmult` step of 0.05 instead of the
  documented 0.5, producing ~10x more combinations than intended, by
  [@sjevelazco](https://github.com/sjevelazco)
- Test suite and vignette runtimes were substantially reduced in
  response to CRAN’s request (test suite ~17 min -\> ~3 min; vignette
  rebuilding ~14 min -\> ~2 min, as measured locally) by shrinking
  hyperparameter grids, dataset sizes, and replicate/fold counts in
  tests and vignettes, and by gating the most expensive/redundant test
  cases behind `skip_on_cran()`, by
  [@sjevelazco](https://github.com/sjevelazco)
- Further reduced test and vignette runtimes after win-builder’s
  automated “Overall checktime \> 10 min” rejection, by gating most
  secondary/edge-case assertions in the `tune_*` and `sdm_uncertainty`
  test files behind `skip_on_cran()` (kept a representative case running
  on CRAN for each), and by shrinking dataset sizes and search-grid
  sizes (`part_sband`, `part_sblock`, `occfilt_env`, `occfilt_geo`)
  further in the vignettes, by
  [@sjevelazco](https://github.com/sjevelazco)
- Test and vignette runtimes were cut again after CRAN’s incoming check
  (tests ~6 min, vignettes ~6 min on win-builder): all but one
  representative `test_that()` block per test file (or none for the
  slowest files) now run only when `NOT_CRAN=true`, and vignettes `v01`,
  `v02`, `v03`, `v05` and `v06` are now precomputed (source in
  `vignettes/*.Rmd.orig`, excluded from the build), by
  [@sjevelazco](https://github.com/sjevelazco)
- Shortened the `\donttest{}` examples of the slowest functions
  (`sdm_varimp`, `fit_max`, `fit_gau`, `tune_*`, `esm_gam`, `fit_gam`,
  `occfilt_geo`, `sdm_summarize`, `esm_max`) by using smaller
  hyperparameter grids, fewer folds/replicates and subsampled datasets
  (~10 min -\> ~3.5 min locally), by
  [@sjevelazco](https://github.com/sjevelazco)
- `fit_raf`: the `ntree` argument was ignored in the no-partition fit
  and in the cross-validation fits (hardcoded to 500); it is now used
  everywhere, reported in
  [\#474](https://github.com/sjevelazco/flexsdm/issues/474), by
  [@sjevelazco](https://github.com/sjevelazco)
- `tune_raf`: the final model now uses the tuned `ntree` instead of the
  `randomForest` default,
  [\#474](https://github.com/sjevelazco/flexsdm/issues/474), by
  [@sjevelazco](https://github.com/sjevelazco)
- `esm_gbm`: `n_trees` and `shrinkage` were hardcoded and are now
  forwarded to `fit_gbm`; `esm_net`: `size` and `decay` are now
  forwarded to `fit_net`,
  [\#474](https://github.com/sjevelazco/flexsdm/issues/474), by
  [@sjevelazco](https://github.com/sjevelazco)
- `esm_svm`: `performance_part` was labelled `esm_gau`; now `esm_svm`.
  `fit_svm` and `fit_max`: partition ids in `performance_part` were
  wrong (earlier folds were relabelled),
  [\#474](https://github.com/sjevelazco/flexsdm/issues/474), by
  [@sjevelazco](https://github.com/sjevelazco)
- Removed the DeepWiki badge from `README.md`, whose URL returned HTTP
  429 in CRAN’s URL check, by
  [@sjevelazco](https://github.com/sjevelazco)
- `sdm_eval`: fixed a bug in the KAPPA metric calculation, by
  [@sjevelazco](https://github.com/sjevelazco)
- `sample_background` and `sample_pseudoabs`: fixed three bugs in
  categorical `maskval`/k-means handling reported in
  [\#472](https://github.com/sjevelazco/flexsdm/issues/472): `maskval`
  on a factor raster was resolved to its row position in the levels
  table instead of its actual category ID; the k-means branch of
  `sample_pseudoabs` masked and overwrote the environmental raster with
  the region raster instead of the other way around; and the internal
  `kf()` helper could fail with “differing number of rows” when the
  environmental raster had no NA cells, by
  [@sjevelazco](https://github.com/sjevelazco)

## flexsdm 1.4.0

- `map_env_dist`: C++ function to calculate Euclidean and Mahalanobis
  distance was implemented, by
  [@sjevelazco](https://github.com/sjevelazco)
  [463](https://github.com/sjevelazco/flexsdm/pull/463)
- `extra_eval`: C++ function to calculate Euclidean and Mahalanobis
  distance was implemented, by
  [@sjevelazco](https://github.com/sjevelazco)
  [463](https://github.com/sjevelazco/flexsdm/pull/463)
- `fit_dom`: C++ function to calculate gower distance was optimized, by
  [@sjevelazco](https://github.com/sjevelazco)
  [463](https://github.com/sjevelazco/flexsdm/pull/463)
- `map_env_dist`: C++ function to calculate gower distance was
  optimized, by [@sjevelazco](https://github.com/sjevelazco)
  [463](https://github.com/sjevelazco/flexsdm/pull/463)
- `msdm_posteriori`: a new argument that assumes unsuitable cells with
  presences as patches. This is only used in the ‘lq’ approach, by
  [@sjevelazco](https://github.com/sjevelazco)
  [461](https://github.com/sjevelazco/flexsdm/pull/461)
- `sdm_eval`: Boyce metric was readapted based on enmSdmX package, by
  [@sjevelazco](https://github.com/sjevelazco)
  [457](https://github.com/sjevelazco/flexsdm/pull/457)
- `msdm_posteriori`: now car return semi-binary suitability maps, by
  [@sjevelazco](https://github.com/sjevelazco)
  [455](https://github.com/sjevelazco/flexsdm/pull/455)
- `part_senv`: new argument to use or not presences-absences coordinate
  to perform environmental partition, by
  [@sjevelazco](https://github.com/sjevelazco)
  [445](https://github.com/sjevelazco/flexsdm/pull/453)

## flexsdm 1.3.9

- `map_env_dist` was improved. Now it is possible to calculate Euclidean
  and Mahalanobis distances in addition to the Gower distance (Domain
  algorithm), by [@sjevelazco](https://github.com/sjevelazco).
- [`sdm_eval()`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.md)
  now calculates three additional performance metrics: the Matthews
  Correlation Coefficient (MCC). the Continuous Ranked Probability Score
  (CRPS), and a prevalence-weighted average of TPR and TNR (W_TPR_TNR),
  by [@sjevelazco](https://github.com/sjevelazco)
  [445](https://github.com/sjevelazco/flexsdm/pull/445).
- [`fit_gbm()`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md)
  and
  [`tune_gbm()`](https://sjevelazco.github.io/flexsdm/reference/tune_gbm.md)
  were improved to prevent the `argument is of length zero` error. This
  was achieved by adding a safeguard that dynamically adjusts
  `n.minobsinnode` based on the size of the training data partition,
  ensuring it remains valid for the `gbm` algorithm, by
  [@sjevelazco](https://github.com/sjevelazco)
  [445](https://github.com/sjevelazco/flexsdm/pull/445).
- `sdm_uncertainty` This new function calculates species distribution
  model uncertainty using a bootstrap procedure, by
  [@sjevelazco](https://github.com/sjevelazco)
  [\#442](https://github.com/sjevelazco/flexsdm/pull/442)
- fit\_ functions now can handle partition = NULL, by
  [@sjevelazco](https://github.com/sjevelazco)
  [\#442](https://github.com/sjevelazco/flexsdm/pull/442)
- A new vignette demonstrating how to model several species and periods
  of time was added by [@mrose048](https://github.com/mrose048)
  [@sjevelazco](https://github.com/sjevelazco)
  [@iminoli-dev](https://github.com/iminoli-dev)
  [438](https://github.com/sjevelazco/flexsdm/pull/438)

## flexsdm 1.3.8

- [`fit_`](https://sjevelazco.github.io/flexsdm/index.html#fit_-functions-family),
  [tune\_](https://sjevelazco.github.io/flexsdm/index.html#tune_-functions-family),
  [esm\_](https://sjevelazco.github.io/flexsdm/index.html#esm_-functions-family)
  and `fit_ensemble` functions now return a performance table for each
  partition and replicate (performance_part), by
  [@sjevelazco](https://github.com/sjevelazco)
  [\#432](https://github.com/sjevelazco/flexsdm/pull/432) and
  [434](https://github.com/sjevelazco/flexsdm/pull/434)
- `p_pdp` was improved to depict exactly training range values when
  projection data are used, by
  [@sjevelazco](https://github.com/sjevelazco)
  [\#429](https://github.com/sjevelazco/flexsdm/pull/429)
- `sample_pseudoabs` a sample approach was implemented based on
  environmental K-means, by [@sjevelazco](https://github.com/sjevelazco)
  [\#410](https://github.com/sjevelazco/flexsdm/pull/410)
- `sample_pseudoabs` K-means step was improved for the
  geoenv_const_kmeans approach, by
  [@sjevelazco](https://github.com/sjevelazco)
  [\#410](https://github.com/sjevelazco/flexsdm/pull/410)
- `p_extra` Now can handle binary maps to plot suitability or
  extrapolation values in the geographical or environmental space, by
  [@sjevelazco](https://github.com/sjevelazco)
  [418](https://github.com/sjevelazco/flexsdm/pull/418)
- `map_env_dist` This new function calculates environmental distance
  between presences and projection data. Only the nearest Gower distance
  was implemented (Domain algorithm), by
  [@sjevelazco](https://github.com/sjevelazco)
  [419](https://github.com/sjevelazco/flexsdm/pull/419/commits/8345526cc87e50d3194030f0ef24f9202bc43a7d)
- `fit_dom` This is a new function to fit and validate the Domain
  algorithm, by [@sjevelazco](https://github.com/sjevelazco)

## flexsdm 1.3.6

- `occfilt_geo` adapted to test different values for the three methods
  by [@sjevelazco](https://github.com/sjevelazco) in
  [386](https://github.com/sjevelazco/flexsdm/pull/386)
- `occfilt_env` now can filter with several bins by
  [@sjevelazco](https://github.com/sjevelazco) in
  [388](https://github.com/sjevelazco/flexsdm/pull/388)
- `occfilt_select` was created and a test protocol was written by
  [@sjevelazco](https://github.com/sjevelazco) in
  [389](https://github.com/sjevelazco/flexsdm/pull/389)
- `occfit_select` included in vignette and website by
  [@sjevelazco](https://github.com/sjevelazco) in
  [390](https://github.com/sjevelazco/flexsdm/pull/390)
- mean change in `sdm_predict` to handle 0 weight value for weighted
  mean by [@sjevelazco](https://github.com/sjevelazco) in
  [391](https://github.com/sjevelazco/flexsdm/pull/391)
- documentation of `occfilt_select` was improved by
  [@sjevelazco](https://github.com/sjevelazco) in
  [392](https://github.com/sjevelazco/flexsdm/pull/392)
- It is possible to tune ntree hyperparameter of random forest algorithm
  by [@sjevelazco](https://github.com/sjevelazco) in
  [393](https://github.com/sjevelazco/flexsdm/pull/393)
- `correct_colinvar` can be used with species points by
  [@sjevelazco](https://github.com/sjevelazco) in
  [394](https://github.com/sjevelazco/flexsdm/pull/394)
- Update pkg_citation.md by [@mrose048](https://github.com/mrose048) in
  [396](https://github.com/sjevelazco/flexsdm/pull/396)
- `sdm_varimp` was added and tested by
  [@mrose048](https://github.com/mrose048) and
  [@sjevelazco](https://github.com/sjevelazco) in
  [399](https://github.com/sjevelazco/flexsdm/pull/399)
- min change in warning message `calib_area` by
  [@sjevelazco](https://github.com/sjevelazco) in
  [400](https://github.com/sjevelazco/flexsdm/pull/400)

## flexsdm 1.3.4-5

- min change in `interp` function by
  [@sjevelazco](https://github.com/sjevelazco) in
  [362](https://github.com/sjevelazco/flexsdm/pull/362)
- min change in `occfilt_geo` by
  [@sjevelazco](https://github.com/sjevelazco) in
  [364](https://github.com/sjevelazco/flexsdm/pull/364)
- Rlof was removed from dependencies by
  [@sjevelazco](https://github.com/sjevelazco) in
  [365](https://github.com/sjevelazco/flexsdm/pull/365)
- min change in `occfilt_geo` by
  [@sjevelazco](https://github.com/sjevelazco) in
  [370](https://github.com/sjevelazco/flexsdm/pull/370)
- new argument in `occfilt_geo` by
  [@sjevelazco](https://github.com/sjevelazco) in
  [371](https://github.com/sjevelazco/flexsdm/pull/371)
- `esm_` functions were improved by
  [@sjevelazco](https://github.com/sjevelazco) in
  [373](https://github.com/sjevelazco/flexsdm/pull/372) &
  [373](https://github.com/sjevelazco/flexsdm/pull/373)
- geographically constraint cell of env_var argument in
  `correct_colinvar` by [@sjevelazco](https://github.com/sjevelazco) in
  [374](https://github.com/sjevelazco/flexsdm/pull/374)
- new argument in `correct_colinvar` by
  [@sjevelazco](https://github.com/sjevelazco) in
  [375](https://github.com/sjevelazco/flexsdm/pull/375)
- `calib_area` was sped up by
  [@sjevelazco](https://github.com/sjevelazco) in
  [377](https://github.com/sjevelazco/flexsdm/pull/377)
- `correct_colivar` was fixed and improved, FA method by
  [@sjevelazco](https://github.com/sjevelazco) in
  [383](https://github.com/sjevelazco/flexsdm/pull/383)
- `occfilt_geo` ‘defined’ method was fixed by
  [@sjevelazco](https://github.com/sjevelazco) in
  [384](https://github.com/sjevelazco/flexsdm/pull/384)

## flexsdm 1.3.4

- it is possible to restrict the cell used to perform collinearity
  reduction analysis to a geographical area smaller than the full extent
  of environmental variables in
  [`correct_clinvar()`](https://sjevelazco.github.io/flexsdm/reference/correct_colinvar.html)
- esm\_ family function was improved and debugged
- `occfilt_geo` has a new argument “rep” to control number o repetition
  to filter occurrences

## flexsdm 1.3.4

- rgeos was removed from dependencies
  [\#356](https://github.com/sjevelazco/flexsdm/pull/356)
- New vignette about how use different tools to explore model
  extrapolation and truncate models was added
  [\#352](https://github.com/sjevelazco/flexsdm/pull/352)
- Univariate and combinatorial extrapolation metric added to
  `extra_eval`.
- Minor bugs were fixed to project PCA for other time periods
  [\#351](https://github.com/sjevelazco/flexsdm/commit/301e241b150d75da4aa01accb3127331ca3bdcb4)
- Best grid raster names was changed to .part in `part_sblock` and
  `part_sband`
- Improvements in `correct_colinvar` to speed up function when using
  maxcell argument
- Improvements in `correct_colinvar` to project PCA for other time
  periods

## flexsdm 1.3.3

- Improvements in `correct_colinvar` It is now possible to sample
  rasters to reduce machine memory and speed up the process
- Improvements in `sdm_predict` It is possible predict model in chunks
  to reduce machine memory
- `p_extra`, `p_pdp`, and `p_bpdp` were fixed
- New function `p_bpdp` Bivariate Partial Dependent Plot
- New function `data_bpdp` Calculate data to construct bivariate partial
  dependence plots
- Improvements in `p_dpd` Calculate data to construct partial dependence
  plots

## flexsdm 1.3.2

- New function `p_extra` Graphical exploration of extrapolation or
  suitability pattern in the environmental and geographical space
- New function `p_pdp` Partial Dependent Plot
- New function `data_pdp` Calculate data to construct partial dependence
  plots

## flexsdm 1.3.1

- New argument “crs” added to function `msdm_posteriori`
- New argument “sp_name” in `sample_background` and `sample_pseudoabs`
- raster, flexclust, ape, and sp were removed from dependencies\
- Functions using CRS data have improved codes
- It is possible use a numeric value to specify threshold in
  `msdm_posteriori`
- `extra_eval` can use tibble or SpatRaster object in env_calib argument
- `extra_truncate` has a new argument to define values used for model
  truncation and documentation was improved. \#
