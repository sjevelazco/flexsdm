# Conventional data partitioning methods

This function provides different conventional (randomized, non-spatial)
partitioning methods based on cross validation folds (kfold, rep_kfold,
and loocv), as well as bootstrap (boot)

## Usage

``` r
part_random(data, pr_ab, method = NULL)
```

## Arguments

- data:

  data.frame. Database with presences, presence-absence, or
  pseudo-absence, records for a given species

- pr_ab:

  character. Column name of "data" with presences, presence-absence, or
  pseudo-absence. Presences must be represented by 1 and absences by 0

- method:

  character. Vector with data partitioning method to be used. Usage
  part=c(method= 'kfold', folds='5'). Methods include:

  - kfold: Random partitioning into k-folds for cross-validation.
    'folds' refers to the number of folds for data partitioning, it
    assumes value \>=1. Usage method = c(method = "kfold", folds = 10).

  - rep_kfold: Random partitioning into repeated k-folds for
    cross-validation. Usage method = c(method = "rep_kfold", folds = 10,
    replicates=10). 'folds' refers to the number of folds for data
    partitioning, it assumes value \>=1. 'replicate' refers to the
    number of replicates, it assumes a value \>=1.

  - loocv: Leave-one-out cross-validation (a.k.a. Jackknife). It is a
    special case of k-fold cross validation where the number of
    partitions is equal to the number of records. Usage method =
    c(method = "loocv").

  - boot: Random bootstrap partitioning. Usage method=c(method='boot',
    replicates='2', proportion='0.7'). 'replicate' refers to the number
    of replicates, it assumes a value \>=1. 'proportion' refers to the
    proportion of occurrences used for model fitting, it assumes a value
    \>0 and \<=1. In this example proportion='0.7' mean that 70% of data
    will be used for model training, while 30% will be used for model
    testing.

## Value

A tibble object with information used in the 'data' argument and
additional columns named .part containing the partition groups. The
rep_kfold and boot method will return as many ".part" columns as
replicated defined. For the rest of the methods, a single .part column
is returned. For kfold, rep_kfold, and loocv partition methods, groups
are defined by integers. In contrast, for boot method, the partition
groups are defined by the characters 'train' and 'test'.

## References

- Fielding, A. H., & Bell, J. F. (1997). A review of methods for the
  assessment of prediction errors in conservation presence/absence
  models. Environmental Conservation, 24(1), 38-49.
  https://doi.org/10.1017/S0376892997000088

## See also

[`part_sblock`](https://sjevelazco.github.io/flexsdm/reference/part_sblock.md),
[`part_senv`](https://sjevelazco.github.io/flexsdm/reference/part_senv.md),
[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md),
[`sample_background`](https://sjevelazco.github.io/flexsdm/reference/sample_background.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("abies")
abies$partition <- NULL
abies <- tibble(abies)

# K-fold method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 10)
)
abies2

# Repeated K-fold method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 10, replicates = 10)
)
abies2

# Leave-one-out cross-validation (loocv) method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "loocv")
)
abies2

# Bootstrap method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 50, proportion = 0.7)
)
abies2
abies2$.part1 %>% table() # Note that for this method .partX columns have train and test words.
} # }
```
