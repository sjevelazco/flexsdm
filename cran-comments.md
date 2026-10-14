## Submission

This is a new submission to CRAN (flexsdm is not currently on CRAN).

## Resubmission

This is a resubmission. In response to the reviewer's comments:

* The title no longer starts with "Tools for" (now "Species Distribution Modeling
  Workflow").
* References describing the methods implemented in the package were added to the
  `Description` field in the requested `authors (year) <doi:...>` form.
* Examples: every exported function now has a small example that is executed
  (no longer wrapped in `\donttest{}`) and runs in well under 5 seconds; the
  longer examples were kept inside `\donttest{}` after the small ones.
* Console output: `print()` / `cat()` calls were replaced by `message()` or
  `warning()` (R/auxiliary_graf.R, R/fit_max.R, R/tune_max.R).
* Writing in the user's file space: `sdm_directory()` no longer has a default
  path (`main_dir` must be provided, previously it defaulted to a sub-folder of
  `getwd()`); the default folder for PCA projections in `correct_colinvar()` is
  now `tempdir()`. Examples, tests and vignettes only write to `tempdir()`.
* Options / graphical parameters: `occfilt_geo()` now restores `options(warn)`
  with `on.exit()`; examples that call `par()` restore it.
* Seeds: functions no longer call `set.seed()` with a fixed number. The seeds used
  by parallel workers in `tune_*()`, `sdm_varimp()` and `sdm_uncertainty()` are
  drawn from the user's random-number stream, so results can be reproduced by
  calling `set.seed()` before running them.
* Functions that used a cluster (`tune_*()`, `sdm_varimp()`, `sdm_uncertainty()`)
  now run sequentially, without starting a cluster, when `n_cores = 1`.

The previous resubmission addressed the incoming check "Overall checktime
> 10 min": most model fitting/tuning tests are skipped on CRAN with
`skip_on_cran()` (they still run when `NOT_CRAN=true`), vignettes `v01`, `v02`,
`v03`, `v05` and `v06` are precomputed (sources in `vignettes/*.Rmd.orig`,
excluded from the build), and the longest examples were shortened.

## Test environments

* local macOS (Apple silicon), R 4.6, `R CMD check --as-cran`
* local Windows 11, R 4.5.1 (`R CMD check --as-cran`)
* win-builder (R-devel): Status 1 NOTE (new submission) in the previous submission

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.
