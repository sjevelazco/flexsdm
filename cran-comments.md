## Submission

This is a new submission to CRAN (flexsdm is not currently on CRAN).

This is a resubmission after the automatic incoming check reported
"Overall checktime 18 min > 10 min" (tests ~357 s, vignette re-building ~369 s).
In response:

* Most `test_that()` blocks, in particular all model fitting/tuning tests
  (`fit_*`, `tune_*`, `esm_*`, `sdm_varimp`, ...), are now skipped on CRAN with
  `skip_on_cran()`; they still run when `NOT_CRAN=true` (e.g. on GitHub Actions).
* Vignettes `v01`, `v02`, `v03`, `v05` and `v06` are now precomputed (sources in
  `vignettes/*.Rmd.orig`, excluded from the build), so only light code is run
  when rebuilding vignettes.
* The `\donttest{}` examples of the slowest functions were shortened (smaller
  grids, fewer folds/replicates, subsampled data).
* The README badge whose URL returned HTTP 429 was removed.

## Test environments

* local Windows 11, R 4.5.1 (`R CMD check --as-cran`)
* win-builder (R-devel): Status 1 NOTE (new submission); tests ~231 s and
  vignette re-building ~60 s before the latest further reduction of the test suite

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.
